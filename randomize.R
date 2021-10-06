kriging_search=function(variable,coordinates,range_search,kappa=0.5,distance= "rdist.earth"){
  ##find range, structure
  ##variable is variable to be kriged (will be normalized).
  ##range_search is set of ranges to try
  ##try again if chosen range is at either end of this range
  ## kappa is smoothness: default exponential 0.5. If ranges too low increase to 1 or 1.5
  ##distance is metric: great circle (default) or euclidean distance= "rdist"
  hold_search=data.frame(range_search,lambda=NA,loglik=NA,eff_df=NA)
  for (j in 1:length(range_search)){
    fit_search = Krig(  x = coordinates,
                        Y = scale(as.vector(variable)),
                        Distance = distance,
                        Covariance = "Matern", 
                        smoothness = kappa,
                        aRange = range_search[j]/sqrt(8*kappa),
                        give.warnings = FALSE
    )
    
    hold_search[j,2:3]=fit_search$lambda.est[1,c(1,5)]
    hold_search[j,4]=fit_search$eff.df
  }
  
  cov_par=hold_search %>% arrange(loglik) %>% slice(1) %>% as.numeric()
  
  if(cov_par[1]==max(range_search)) stop("Your range search values are too low: try higher ones or increase smoothness kappa.")
  if(cov_par[1]==min(range_search)) stop("Your range search values are too high: try lower ones or reduce smoothness kappa.")
  
  Effective_Range=1.6*sqrt(8*kappa)*cov_par[1]  #in kilometres
  
  results=data.frame(Range=cov_par[1],Structure=1/(1+cov_par[2]),Effective_Range=Effective_Range,
                     Noise_to_Signal=cov_par[2],Likelihood=cov_par[3],df=cov_par[4]/nrow(coordinates))
  
  return(results)
}


#########randomize x

x_randomize=function(study,frm,Range_Search,poly_degree=2,kappa=0.5,
                     nSim=1000 ,
                     lower_cutoff= -100, upper_cutoff= 100){
  #range search is effective range in km
  
  ###################first turn variables into explan and dependent
  frm=str_remove_all(frm," ")
  frm=str_remove_all(frm,"\\\n")
  frm_spl=str_split(frm,"~")
  dep_var=frm_spl[[1]][1]
  
  if (str_detect(frm,"\\+")) {
    ex=str_split(frm_spl[[1]][2],"\\+",n=2)
    explan_var=ex[[1]][1]
    rhs_vars=paste("+",ex[[1]][2])
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
  nearest = nb2listw(nearest, style = "W")
  Moran_z = as.vector(lm.morantest(ols, listw = nearest)$statistic) 
  
  
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

  Output_x=data.frame(
    p_exact_x=p_exact,
    p_orig,
    t_orig,
    direction_R2=direction_R2,
    Effective_Range=res_params$Effective_Range,
    Structure=Structure,
    kappa=kappa,
    Moran=Moran_z,
    Likelihood=res_params$Likelihood
  )
  
  return(list(Output_x=Output_x,
              Residuals=Residuals,
              t_noise_x=sim_res$t_noise))
}

y_randomize=function(study,frm,Range_Search,
                     poly_degree=2,kappa=0.5,
                     nSim=1000 ,
                     lower_cutoff= -100, upper_cutoff= 100){
  
  ###################first turn variables into explan and dependent
  frm=str_remove_all(frm," ")
  frm=str_remove_all(frm,"\\\n")
  frm_spl=str_split(frm,"~")
  dep_var=frm_spl[[1]][1]
  
  if (str_detect(frm,"\\+")) {
    ex=str_split(frm_spl[[1]][2],"\\+",n=2)
    explan_var=ex[[1]][1]
    rhs_vars=paste("+",ex[[1]][2])
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
  nearest = nb2listw(nearest, style = "W")
  Moran_z = as.vector(lm.morantest(ols, listw = nearest)$statistic) 
  
  
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
  

  Output_y=data.frame(
    p_exact_y=p_exact,
    p_orig,
    t_orig,
    direction_R2=direction_R2,
    Effective_Range=res_params$Effective_Range,
    Structure=Structure,
    kappa=kappa,
    Moran=Moran_z,
    Likelihood=res_params$Likelihood
  )
  
  return(list(Output_y=Output_y,
              Residuals=Residuals,
              t_noise_y=sim_res$t_noise)
         )
}

############Average correlation for Mueller-West estimation
avg_correlation=function(study,frm,kappa=0.5,Range_Search,
                         lower_cutoff= -100, upper_cutoff= 100){
  #range search is effective range in km
  
  ###################first turn variables into explan and dependent
  frm=str_remove_all(frm," ")
  frm=str_remove_all(frm,"\\\n")
  frm_spl=str_split(frm,"~")
  dep_var=frm_spl[[1]][1]
  
  if (str_detect(frm,"\\+")) {
    ex=str_split(frm_spl[[1]][2],"\\+",n=2)
    explan_var=ex[[1]][1]
    rhs_vars=paste("+",ex[[1]][2])
  } else{
    explan_var=frm_spl[[1]][2]
    rhs_vars=""
  }
  
  study=study %>% mutate(dep_var=study[,dep_var],explan_var=study[,explan_var])
  
  ################################################################
  study=study%>%
    mutate(abs_Y=abs(Y),X1=X*pi/180,abs_Y1=abs_Y*pi/180)  #radians
  
  fm_ols=as.formula(paste("dep_var~explan_var",rhs_vars,sep=""))
  
  
  
  ##original ols formulae 
  Residuals=lm(fm_ols,study)$residuals
  
  Residuals[Residuals<lower_cutoff]=lower_cutoff
  Residuals[Residuals>upper_cutoff]=upper_cutoff
  
  #############Correlation parameters
  resids=scale(Residuals)
  Coords=study %>% dplyr::select(X,Y) %>% as.matrix()
  Range_Search=Range_Search/(1.6*sqrt(8*kappa))   #convert to kilometres and turn effective range into range
  res_params=kriging_search(resids,Coords,range_search=Range_Search,kappa=kappa)
  Range=as.numeric(res_params$Range)
  Structure=as.numeric(res_params$Structure)  #this is weight rho between systematic correlation and spatial noise
  
  
  KL=Structure*fields::Matern(rdist.earth(x1=Coords),
                              range=Range,smoothness=kappa)+
    diag(nrow(Coords))*(1-Structure)
  diag(KL)=NA
  avg_corr=mean(KL,na.rm=T)
  
  return(avg_corr)
}
# aa=x_randomize(study,rhs_vars,Range_Search=seq(100,1000,by=100),kappa=0.5,nSim=1000)
# aa$Output
# # densityPlot(aa$t_noise)
# # densityPlot(aa$Residuals)
# # 
# # 
# study %>%
#   mutate(explan=santoku::chop_deciles(explan_var)) %>%
#   ggplot(aes(X,Y,colour=explan))+geom_point(shape=15,size=1.75)+
#   theme_minimal(base_size = 16) +
#   scale_colour_brewer(palette = "RdYlBu",direction= -1)+
#   theme(legend.position = "none")+
# labs(subtitle="Explanatory Variable")
# 
# study %>%
#   mutate(dep=santoku::chop_deciles(dep_var)) %>%
#   ggplot(aes(X,Y,colour=dep))+geom_point(shape=15,size=1.75)+
#   scale_colour_brewer(palette = "RdYlBu")+
#   theme_minimal(base_size = 16) +
#   theme(legend.position = "none")+
#   labs(subtitle="Dependent Variable")
# 
# study %>% ggplot(aes(explan_var,dep_var))+
#   geom_point(size=2)+
#   theme_bw()
# 
# ggplot(study, aes(explan_var, dep_var, colour = factor(region))) + geom_point(size = 2) +
#   scale_color_viridis_d(direction = -1) + theme_minimal(base_size = 15) +
#   theme(legend.position = "none") + labs(x = "Mortality", y = "Property Rights",
#                                          title = "Mortality and Property Rights", subtitle = "Shaded by WB region.")
