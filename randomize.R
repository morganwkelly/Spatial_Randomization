kriging_search=function(variable,coordinates,range_search,kappa=0.5,distance= "rdist.earth"){
  ##find range, structure
  ##variable is variable to be kriged (will be normalized).
  ##range_search is set of ranges to try
  ##try again if chosen range is at either end of this range
  ## kappa is smoothness: default exponential 0.5. If ranges too low increase to 1 or 1.5
  ##distance is metric: great circle (default) or euclidean distance= "rdist"
  hold_search=data.frame(range_search,lambda=NA,loglik=NA,converge=NA,eff_df=NA)
  for (j in 1:length(range_search)){
    fit_search = Krig(  x = coordinates,
                        Y = scale(as.vector(variable)),
                        Distance = distance,
                        Covariance = "Matern", 
                        smoothness = kappa,
                        aRange = range_search[j]/sqrt(8*kappa),
                        give.warnings = F
    )
    
    hold_search[j,2:4]=fit_search$lambda.est[6,c(1,5,6)]
    hold_search[j,5]=fit_search$eff.df
  }
  
  cov_par=hold_search %>% arrange(loglik) %>% slice(1) %>% as.numeric()
  
   if(cov_par[1]==max(range_search)) stop("Your range search values are too low: try higher ones or increase smoothness kappa.")
   if(cov_par[1]==min(range_search)) stop("Your range search values are too high: try lower ones or reduce smoothness kappa.")
  # 
  Effective_Range=1.6*sqrt(8*kappa)*cov_par[1]  #in kilometres
  
  results=data.frame(Range=cov_par[1],Structure=1/(1+cov_par[2]),Effective_Range=Effective_Range,
                     Noise_to_Signal=cov_par[2],Likelihood=cov_par[3],
                     df=cov_par[5]/nrow(coordinates))
  
  return(results)
}


#########randomize x

x_randomize=function(study,frm,var_num=1,Range_Search,poly_degree=2,kappa=0.5,
                     nSim=2000 ,
                     lower_cutoff= -100, upper_cutoff= 100){
  #range search is effective range in km
  
  ###################first turn variables into explan and dependent
  frm=str_remove_all(frm," ")
  frm=str_remove_all(frm,"\\\n")
  frm_spl=str_split(frm,"~")
  dep_var=frm_spl[[1]][1]
  
  if (str_detect(frm,"\\+")) {
    ex=str_split(frm_spl[[1]][2],"\\+")[[1]]
    explan_var=ex[var_num]
    rhs_vars=paste("+",paste(ex[-var_num],collapse="+"))
  } else{
    explan_var=frm_spl[[1]][2]
    rhs_vars=""
  }
  study=as.data.frame(study)
  study=study %>% mutate(dep_var=study[,dep_var],explan_var=study[,explan_var])
  

  
  
  ################################################################
  study=study%>%
    mutate(abs_Y=abs(Y),X1=X*pi/180,abs_Y1=abs_Y*pi/180)  #radians

  fm_ols=as.formula(paste("dep_var~explan_var",rhs_vars,sep=""))
  fm_xsim=as.formula(paste("dep_var~sim_x",rhs_vars,sep=""))
  fm_orthog=as.formula(paste("explan_var~",rhs_vars,sep=""))

  
 ##original ols formulae 
  lm_rob=lm_robust(fm_ols,study,se_type="stata")
  rob_sum=summary(lm_rob)
  #rob_sum
  coef_orig=rob_sum$coef[2,1]
  se_orig=rob_sum$coef[2,2]
  t_orig=rob_sum$coef[2,3]
  p_orig=rob_sum$coef[2,4]
  R2=rob_sum$r.squared
  
  ############orthog variable
  lm_orthog=lm(fm_orthog,data=study)
  explan_orthog=lm_orthog$residuals
  
  lm_res=lm(explan_orthog~poly(abs_Y1,poly_degree)+poly(X1,poly_degree),data=study)
  summary(lm_res)
  Residuals=lm_res$residuals
  direction_R2=summary(lm_res)$r.squared  #explan power of spatial direction
  Residuals[Residuals<lower_cutoff]=lower_cutoff
  Residuals[Residuals>upper_cutoff]=upper_cutoff

  ###########compute moran stat for spatial autocorrelation
  Coords=study %>% dplyr::select(X,Y) %>% as.matrix()
  ols=lm(fm_ols,study)
  study_spatial = SpatialPoints(coords = Coords)
  proj4string(study_spatial) = CRS("+proj=longlat +datum=WGS84")
  
  nearest = knn2nb(knearneigh(study_spatial, k = 5, longlat = T))   #k nearest neighbours
  Moran_p = lm.morantest(ols, listw=nb2listw(nearest,style="W"))$p.value[1,1] 
  
  
  #############Correlation parameters
  Range_Search=Range_Search/(1.6*sqrt(8*kappa))   #convert to kilometres and turn effective range into range
  res_params=kriging_search(Residuals,Coords,range_search=Range_Search,kappa=kappa)
   Range=as.numeric(res_params$Range)
  Structure=as.numeric(res_params$Structure)  #this is weight rho between systematic correlation and spatial noise

 
  KL=Structure*fields::Matern(rdist.earth(x1=Coords),
                              range=Range,smoothness=kappa)+
    diag(nrow(Coords))*(1-Structure)
  
  L_res=t(chol(KL))
  
  nTot=nSim
  set.seed(1234)

  xSim=L_res%*%matrix(rnorm(nTot*nrow(Coords)),ncol=nTot)
  xSim=lm_res$fitted+sd(Residuals)*xSim
  
  
  sim_res=data.frame(coef=rep(NA,nSim),t_noise=rep(NA,nSim))
  for (i in 1:nSim){
      sim_x=xSim[,i]
      rob_se=lm_robust(fm_xsim,study,se_type="stata")
      sim_res[i,1:2]=summary(rob_se)$coef[2,c(1,3)]
    }

  p_exact=sum(abs(sim_res$t_noise)>abs(t_orig))/nSim
  
  se_exact=abs(coef_orig/qnorm(p_exact/2))  #implicit SE
 
  Output_x=data.frame(
    explan_var,
    p_exact_x=p_exact,
    p_orig,
    t_orig=abs(t_orig),
    coefficient=coef_orig,
    se_exact_x=se_exact,
    se_orig,
    R2,
    N=nrow(study),
    direction_R2_x=direction_R2,
    Effective_Range_x=res_params$Effective_Range,
    Structure_x=Structure,
    kappa_x=kappa,
    Moran_p,
    Likelihood_x=res_params$Likelihood,
    df_x=res_params$df
  )
  
  return(list(Output_x=Output_x,
              Residuals=Residuals,
              t_noise_x=sim_res$t_noise))
}

y_randomize=function(study,frm,Range_Search,var_num=1,
                     poly_degree=2,kappa=0.5,
                     nSim=1000 ,
                     lower_cutoff= -100, upper_cutoff= 100){
  
  ###################first turn variables into explan and dependent
  frm=str_remove_all(frm," ")
  frm=str_remove_all(frm,"\\\n")
  frm_spl=str_split(frm,"~")
  dep_var=frm_spl[[1]][1]
  
  if (str_detect(frm,"\\+")) {
    ex=str_split(frm_spl[[1]][2],"\\+")[[1]]
    explan_var=ex[var_num]
    rhs_vars=paste("+",paste(ex[-var_num],collapse="+"))
  } else{
    explan_var=frm_spl[[1]][2]
    rhs_vars=""
  }
  study=as.data.frame(study)
  study=study %>% mutate(dep_var=study[,dep_var],explan_var=study[,explan_var])
  
  ################################################################
  
  study=study%>%
    mutate(abs_Y=abs(Y),X1=X*pi/180,abs_Y1=abs_Y*pi/180)    #radians
  
  fm_ols=as.formula(paste("dep_var~explan_var",rhs_vars,sep=""))
  fm_ols
  fm_ysim=as.formula(paste("sim_y~explan_var",rhs_vars,sep=""))
  fm_yorthog=as.formula(paste("dep_var~",rhs_vars,sep=""))
  
  ##original ols formulae 
  lm_rob=lm_robust(fm_ols,study,se_type="stata")
  rob_sum=summary(lm_rob)
  #rob_sum
  coef_orig=rob_sum$coef[2,1]
  se_orig=rob_sum$coef[2,2]
  t_orig=rob_sum$coef[2,3]
  p_orig=rob_sum$coef[2,4]
  R2=rob_sum$r.squared
  
  ############orthog variable
  lm_orthog=lm(fm_yorthog,data=study)
  dep_orthog=lm_orthog$residuals
  
  lm_res=lm(dep_orthog~poly(abs_Y1,poly_degree)+poly(X1,poly_degree),data=study)
  summary(lm_res)
  Residuals=lm_res$residuals
  direction_R2=summary(lm_res)$r.squared  #explan power of spatial direction
  
  Residuals[Residuals<lower_cutoff]=lower_cutoff
  Residuals[Residuals>upper_cutoff]=upper_cutoff
  
  
  ###########compute moran stat for spatial autocorrelation
  Coords=study %>% dplyr::select(X,Y) %>% as.matrix()
  ols=lm(fm_ols,study)
  study_spatial = SpatialPoints(coords = Coords)
  proj4string(study_spatial) = CRS("+proj=longlat +datum=WGS84")
  
  nearest = knn2nb(knearneigh(study_spatial, k = 5, longlat = T))   #k nearest neighbours
  Moran_p = lm.morantest(ols, listw=nb2listw(nearest,style="W"))$p.value[1,1] 
  
  
  
  #############Correlation parameters
  Coords=study %>% dplyr::select(X,Y) %>% as.matrix()
  Range_Search=Range_Search/(1.6*sqrt(8*kappa))   #convert effective range in km to range in miles
  res_params=kriging_search(Residuals,Coords,range_search=Range_Search,kappa=kappa)
  Range=as.numeric(res_params$Range)
  Structure=as.numeric(res_params$Structure)  #this is weight rho between systematic correlation and spatial noise
  
  
  KL=Structure*fields::Matern(rdist.earth(x1=Coords),
                              range=Range,smoothness=kappa)+
    diag(nrow(Coords))*(1-Structure)
  
  L_res=t(chol(KL))
  
  nTot=nSim
  set.seed(1234)
  
  ySim=L_res%*%matrix(rnorm(nTot*nrow(Coords)),ncol=nTot)
  ySim=apply(ySim,2,scale)  #scale
  ySim=lm_res$fitted+sd(Residuals)*ySim
  ySim=lm_orthog$fitted.values+ySim
  
  
  sim_res=data.frame(coef=rep(NA,nSim),t_noise=rep(NA,nSim))
  for (i in 1:nSim){
    sim_y=ySim[,i]
    lm_ysim=lm(fm_ysim,data=study)
    rob_se=lm_robust(fm_ysim,study,se_type="stata")
    sim_res[i,1:2]=summary(rob_se)$coef[2,c(1,3)]
  }
  
  
  p_exact=sum(abs(sim_res$t_noise)>abs(t_orig))/nSim
  
  se_exact=abs(coef_orig/qnorm(p_exact/2))  #implicit SE

  Output_y=data.frame(
    explan_var,
    p_exact_y=p_exact,
    p_orig,
    t_orig=abs(t_orig),
    coefficient=coef_orig,
    se_exact_y=se_exact,
    se_orig,
    R2,
    N=nrow(study),
    direction_R2_y=direction_R2,
    Effective_Range_y=res_params$Effective_Range,
    Structure_y=Structure,
    kappa_y=kappa,
    Moran_p,
    Likelihood_y=res_params$Likelihood,
    df_y=res_params$df
  )
  
  return(list(Output_y=Output_y,
              Residuals=Residuals,
              t_noise_y=sim_res$t_noise)
         )
}



randomized_crit_vals=function(study,frm,N_vars,range,sims=1000){
########make table of randomized values  
  out_x=list()
  out_y=list()
  
  for (i in 1:N_vars){
    n_var=i
    x_1=x_randomize(study,frm,var_num=n_var,
                    Range_Search=range,poly_degree=2,kappa=0.5,nSim=sims)
    x_output=x_1$Output_x
    out_x[[i]]=x_output
    
    y_1=y_randomize(study,frm,var_num=n_var,
                    Range_Search=range,poly_degree=2,kappa=0.5,nSim=sims)
    y_output=y_1$Output_y
    #y_output$MW_Corr=NA    #avg_correlation(study,frm,var_num=n_var,Range_Search=range)
    out_y[[i]]=y_output
  }
  
  out_x=plyr::ldply(out_x)
  out_y=plyr::ldply(out_y)
  return(list(x=out_x,y=out_y))
}

make_model=function(eqn,rand_est){
  # Generate model summary output for randomized results for modelsummary
  # eqn is estimated eqn
  # rand_est is return from randomize summary
  a1=str_replace_all(c(eqn),"\\+",",")
  a2=str_split(a1,"\\~")[[1]][2]
  a3=str_split(a2,",")[[1]]
  a3=a3[1:nrow(rand_est$x)]
  
  ti <- data.frame(
    term = a3,
    estimate = rand_est[[1]]$coefficient,
    std.error=rand_est[[1]]$se_orig,
    random.std.error = apply(cbind(rand_est[[1]]$se_exact_x,rand_est[[2]]$se_exact_y,    #include original SE in case exact p=0
                            rand_est[[1]]$se_orig),
                      1,max)
  )
  
  gl <- data.frame(
    N=rand_est[[1]]$N[1],
    R2=rand_est[[1]]$R2[1],
    Moran_p = rand_est[[1]]$Moran[1])
  
  mod <- list(
    tidy = ti,
    glance = gl)
  class(mod) <- "modelsummary_list"
  return(mod)
}


iv_randomize=function(study,frm,Range_Search,var_num=1,
                      poly_degree=2,kappa=0.5,
                      nSim=1000 ,
                      lower_cutoff= -100, upper_cutoff= 100){
  
  ###################first turn variables into explan and dependent
  frm=str_remove_all(frm," ")
  frm=str_remove_all(frm,"\\\n")
  frm_iv=str_split(frm,"\\|")[[1]][2]
  frm=str_split(frm,"\\|")[[1]][1]
  frm_spl=str_split(frm,"~")
  dep_var=frm_spl[[1]][1]
  
  if (str_detect(frm,"\\+")) {
    ex=str_split(frm_spl[[1]][2],"\\+")[[1]]
    explan_var=ex[var_num]
    rhs_vars=paste("+",paste(ex[-var_num],collapse="+"))
  } else{
    explan_var=frm_spl[[1]][2]
    rhs_vars=""
  }
  
  iv=str_split(frm_iv,"\\+")[[1]]
  if(var_num>1){
    iv[var_num]="explan_var"
    iv_vars=paste(iv,collapse="+")
    iv_vars=paste("|",iv_vars,sep="")
  }else{
    iv_vars=paste("|",frm_iv) 
  }
  study=as.data.frame(study)
  study=study %>% mutate(dep_var=study[,dep_var],explan_var=study[,explan_var])
  
  ################################################################
  
  study=study%>%
    mutate(abs_Y=abs(Y),X1=X*pi/180,abs_Y1=abs_Y*pi/180)    #radians
  
  fm_ols=as.formula(paste("dep_var~explan_var",rhs_vars,iv_vars,sep=""))
  fm_ols
  fm_ysim=as.formula(paste("sim_y~explan_var",rhs_vars,iv_vars,sep=""))
  fm_yorthog=as.formula(paste("dep_var~",rhs_vars,iv_vars,sep=""))
  
  ##original ols formulae 
  lm_rob=iv_robust(fm_ols,study,se_type="stata")
  rob_sum=summary(lm_rob)
  rob_sum
  coef_orig=rob_sum$coef[2,1]
  se_orig=rob_sum$coef[2,2]
  t_orig=rob_sum$coef[2,3]
  p_orig=rob_sum$coef[2,4]
  R2=rob_sum$r.squared
  
  ############orthog variable
  lm_orthog=lm(fm_yorthog,data=study)
  dep_orthog=lm_orthog$residuals
  
  lm_res=lm(dep_orthog~poly(abs_Y1,poly_degree)+poly(X1,poly_degree),data=study)
  #summary(lm_res)
  Residuals=lm_res$residuals
  direction_R2=summary(lm_res)$r.squared  #explan power of spatial direction
  
  Residuals[Residuals<lower_cutoff]=lower_cutoff
  Residuals[Residuals>upper_cutoff]=upper_cutoff
  
  ############diagnostics for IV
  iv1=ivreg(fm_ols,data=study)
  iv_diags=summary(iv1,diagnostics = TRUE)$diagnostics
  Weak_Instruments=iv_diags[1,4]
  Wu_Hausman=iv_diags[2,4]
  Sargan=iv_diags[3,4]
  
  ###########compute moran stat for spatial autocorrelation
  Coords=study %>% dplyr::select(X,Y) %>% as.matrix()
  ols=lm(fm_ols,study)
  study_spatial = SpatialPoints(coords = Coords)
  proj4string(study_spatial) = CRS("+proj=longlat +datum=WGS84")
  
  nearest = knn2nb(knearneigh(study_spatial, k = 5, longlat = T))   #k nearest neighbours
  nearest = nb2listw(nearest, style = "W")
  Moran_p = as.vector(moran.test(iv1$residuals,listw=nb2listw(contig,style="W"))$p) 
  
  
  #############Correlation parameters
  Coords=study %>% dplyr::select(X,Y) %>% as.matrix()
  Range_Search=Range_Search/(1.6*sqrt(8*kappa))   #convert effective range in km to range in miles
  res_params=kriging_search(Residuals,Coords,range_search=Range_Search,kappa=kappa)
  Range=as.numeric(res_params$Range)
  Structure=as.numeric(res_params$Structure)  #this is weight rho between systematic correlation and spatial noise
  
  
  KL=Structure*fields::Matern(rdist.earth(x1=Coords),
                              range=Range,smoothness=kappa)+
    diag(nrow(Coords))*(1-Structure)
  
  L_res=t(chol(KL))
  
  nTot=nSim
  set.seed(1234)
  
  ySim=L_res%*%matrix(rnorm(nTot*nrow(Coords)),ncol=nTot)
  ySim=apply(ySim,2,scale)  #scale
  ySim=lm_res$fitted+sd(Residuals)*ySim
  ySim=lm_orthog$fitted.values+ySim
  
  
  sim_res=data.frame(coef=rep(NA,nSim),t_noise=rep(NA,nSim))
  for (i in 1:nSim){
    sim_y=ySim[,i]
    lm_ysim=lm(fm_ysim,data=study)
    rob_se=iv_robust(fm_ysim,study,se_type="stata")
    sim_res[i,1:2]=summary(rob_se)$coef[2,c(1,3)]
  }
  
  
  p_exact=sum(abs(sim_res$t_noise)>abs(t_orig))/nSim
  
  se_exact=abs(coef_orig/qnorm(p_exact/2))  #implicit SE
  
  
 
  
  Output_IV=data.frame(
    explan_var,
    p_exact_y=p_exact,
    p_orig,
    t_orig=abs(t_orig),
    coefficient=coef_orig,
    se_exact_y=se_exact,
    se_orig,
    R2,
    N=iv1$n,
    direction_R2_y=direction_R2,
    Effective_Range_y=res_params$Effective_Range,
    Structure_y=Structure,
    kappa_y=kappa,
    Moran_p,
    Weak_Instruments,
    Wu_Hausman,
    Sargan,
    Likelihood_y=res_params$Likelihood,
    df_y=res_params$df
  )
  
  return(list(Output_IV=Output_IV,
              Residuals=Residuals,
              t_noise_y=sim_res$t_noise)
  )
}

randomized_iv_crit_vals=function(study,frm,N_vars,range,sims=1000){
  ########make table of randomized values  
  out_iv=list()
  
  for (i in 1:N_vars){
    n_var=i
    x_iv=iv_randomize(study,frm,var_num=n_var,
                      Range_Search=range,poly_degree=2,kappa=0.5,nSim=sims)
    x_output=x_iv$Output_IV
    out_iv[[i]]=x_output
  }
  
  out_iv=plyr::ldply(out_iv)
  return(out_iv)
}


###########note: input here is df, not list
make_iv_model=function(eqn,rand_est){
  # Generate model summary output for randomized results for modelsummary
  # eqn is estimated eqn
  # rand_est is return from randomize summary
  a1=str_replace_all(c(eqn),"\\+",",")
  a2=str_split(a1,"\\~")[[1]][2]
  a3=str_split(a2,",")[[1]]
  a3=a3[1:nrow(rand_est)]
  
  ti <- data.frame(
    term = a3,
    estimate = rand_est$coefficient,
    std.error=rand_est$se_orig,
    random.std.error = apply(cbind(rand_est$se_exact_y,    #include original SE in case exact p=0
                                   rand_est$se_orig),
                             1,max)
  )
  
  gl <- data.frame(
    N=rand_est$N[1],
    Moran_p = rand_est$Moran[1],
    Weak_Instruments = rand_est$Weak_Instruments[1],
    Wu_Hausman = rand_est$Wu_Hausman[1],
    Sargan = rand_est$Sargan[1]
  )
  
  mod <- list(
    tidy = ti,
    glance = gl)
  class(mod) <- "modelsummary_list"
  return(mod)
}



