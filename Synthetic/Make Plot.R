  # Make Figures for the Synthetic Version of SC Prenatal Data
  # Import: MCMC Samples (File: Simulation-ZIBBSVC.Rda) 
  # Note: Results can be differed from the manuscript (synthetic data)
  # Date: 10-22-2024
  #------------------------------------------------------------------------------#
  library(aod)
  library(VGAM)       # For rbetabinom()
  library(mcmcplots)
  library(spam)
  library(dplyr)
  library(gamlss)     # For dZIBB()
  library(tidyr)
  library(tmap)
  library(tigris)
  library(RColorBrewer)
  library(ggpubr)
  library(ggplot2)
  #------------------------------------------------------------------------------#
  rm(list=ls())
  #------------------------------------------------------------------------------#
  
  #---------------------#
  # Load Synthetic Data #
  #---------------------#
  dir.sav<-"...\\"
  load(file=paste(dir.sav,"Synthetic-CRF-Data.RData",sep=""))
  dsvi=read.csv("...\\SC_SVI.csv")

  #---------------#
  # Load SVI Data #
  #---------------#
  n<-46 				                   # Number of SC county
  nis<-dsyn%>%group_by(countyid)%>%summarise(nis=n())%>%pull(nis)  # Number of pregnancies in each county
  ntrial<-5                        # Total number of CRFs
  id<-rep(1:n,nis)                 # County ID
  N<-length(id)                    # Number of total pregnancies
  
  # Covariate
  t<-dsyn$timeqr-5                    # Time variable from Q1,2020 to Q4, 2020 (centered)
  nhb<-dsyn$nhb                       # NHB indicator
  int<-t*nhb                       # Interaction (time x race)
  agec<-dsyn$age-mean(dsyn$age)          # Maternal age
  edu1<-dsyn$edu1                     # Education (< high school)
  edu2<-dsyn$edu2                     # high school
  edu3<-dsyn$edu3                     # some college
  medicaid<-dsyn$medicaid             # Medicaid use during pregnancy
  
  # Covaraite -- County
  svi<-dsvi$svi_overall_rank       # SVI, percentile-based
  svis<-scale(svi)                 # standardize
  SVIS<-rep(svis,nis)              
  
  # Design matrix
  X<-cbind(1,t,nhb,int,SVIS,agec,edu1,edu2,edu3,medicaid)
  p<-ncol(X)
  
  # Response
  y<-dsyn$y                           # number of CRFs
  
  #-------------------#
  # Load MCMC Samples #
  #-------------------#
  load(file=paste(dir.sav,"Synthetic-MCMCSampler.Rda",sep=""))           # Load MCMC samples
  Alpha<-samples$Alpha
  Beta<-samples$Beta
  Rho<-samples$Rho
  S2Phi<-samples$S2Phi
  Phi1s<-samples$Phi1s
  Phi2s<-samples$Phi2s
  Phi3s<-samples$Phi3s
  Phi4s<-samples$Phi4s
  Phi5s<-samples$Phi5s
  Phi6s<-samples$Phi6s
  Phi7s<-samples$Phi7s
  Phi8s<-samples$Phi8s
  
  #--------------------#
  # MCMC Sampler Parm. #
  #--------------------#
  nsim<-20000                   # Number of MCMC Iterations
  thin<-20				              # Thinnisng interval
  burn<-nsim/2   	              # Burnisn
  lastit<-(nsim-burn)/thin     	# Last stored value
  
  #------------------------------------------------------------------------------#
  # Odds Ratio #
  #------------#

  #----------------------------------------------------#
  # Marginal CRF (Average County, all random effect=0) #
  #----------------------------------------------------#
  
  # NOTE: Q1,2020 (t=-4), Q4,2020 (t=-1), Q4,2021 (t=2), Q4,2022 (t=7)
  
  # Average county mean No. of CRFs at Q1,2020
  te<--4
  eta1<-Alpha[,1]+te*Alpha[,2]
  mu1_t1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+te*Beta[,2]
  mu2_t1<-exp(eta2)/(1+exp(eta2))
  mean(mu1_t1*mu2_t1*5);quantile(mu1_t1*mu2_t1*5,c(.025,.975)) # Predicted number of CRFs
  
  #---------------------------#
  # Marginal Mean Ratio (MMR) #
  #---------------------------#
  
  ### SVI = 1 sd increase
  mean(exp(Alpha[,5]))
  quantile(exp(Alpha[,5]),c(0.025,0.975))
  
  eta1<-Alpha[,1]
  mu1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]
  mu2<-exp(eta2)/(1+exp(eta2))
  
  eta1<-Alpha[,1]+Alpha[,5]
  mu1_svi<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+Beta[,5]
  mu2_svi<-exp(eta2)/(1+exp(eta2))
  
  mean(mu1_svi*mu2_svi/mu1/mu2);quantile(mu1_svi*mu2_svi/mu1/mu2,c(.025,.975))
  
  ### Age = 1 unit increase 
  eta1<-Alpha[,1]
  mu1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]
  mu2<-exp(eta2)/(1+exp(eta2))
  
  eta1<-Alpha[,1]+Alpha[,6]
  mu1_age<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+Beta[,6]
  mu2_age<-exp(eta2)/(1+exp(eta2))
  
  mean(mu1_age*mu2_age/mu1/mu2);quantile(mu1_age*mu2_age/mu1/mu2,c(.025,.975))
  
  ### < High School
  eta1<-Alpha[,1]
  mu1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]
  mu2<-exp(eta2)/(1+exp(eta2))
  
  eta1<-Alpha[,1]+Alpha[,7]
  mu1_lhs<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+Beta[,7]
  mu2_lhs<-exp(eta2)/(1+exp(eta2))
  
  mean(mu1_lhs*mu2_lhs/mu1/mu2);quantile(mu1_lhs*mu2_lhs/mu1/mu2,c(.025,.975))
  
  ### High School 
  eta1<-Alpha[,1]
  mu1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]
  mu2<-exp(eta2)/(1+exp(eta2))
  
  eta1<-Alpha[,1]+Alpha[,8]
  mu1_hs<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+Beta[,8]
  mu2_hs<-exp(eta2)/(1+exp(eta2))
  
  mean(mu1_hs*mu2_hs/mu1/mu2);quantile(mu1_hs*mu2_hs/mu1/mu2,c(.025,.975))
  
  ### Some College 
  eta1<-Alpha[,1]
  mu1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]
  mu2<-exp(eta2)/(1+exp(eta2))
  
  eta1<-Alpha[,1]+Alpha[,9]
  mu1_sc<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+Beta[,9]
  mu2_sc<-exp(eta2)/(1+exp(eta2))
  
  mean(mu1_sc*mu2_sc/mu1/mu2);quantile(mu1_sc*mu2_sc/mu1/mu2,c(.025,.975))
  
  ### Medicaid 
  eta1<-Alpha[,1]
  mu1<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]
  mu2<-exp(eta2)/(1+exp(eta2))
  
  eta1<-Alpha[,1]+Alpha[,10]
  mu1_med<-exp(eta1)/(1+exp(eta1))
  eta2<-Beta[,1]+Beta[,10]
  mu2_med<-exp(eta2)/(1+exp(eta2))
  
  mean(mu1_med*mu2_med/mu1/mu2);quantile(mu1_med*mu2_med/mu1/mu2,c(.025,.975))
  
  
  #------------------------------------------------------------------------------#
  # Make Map #
  #----------#
  
  mphi1<-colMeans(Phi1s)
  mphi2<-colMeans(Phi2s)
  mphi3<-colMeans(Phi3s)
  mphi4<-colMeans(Phi4s)
  mphi5<-colMeans(Phi5s)
  mphi6<-colMeans(Phi6s)
  mphi7<-colMeans(Phi7s)
  mphi8<-colMeans(Phi8s)
  
  qphi1<-apply(Phi1s,2,quantile,c(.025,.975))
  qphi2<-apply(Phi2s,2,quantile,c(.025,.975))
  qphi3<-apply(Phi3s,2,quantile,c(.025,.975))
  qphi4<-apply(Phi4s,2,quantile,c(.025,.975))
  qphi5<-apply(Phi5s,2,quantile,c(.025,.975))
  qphi6<-apply(Phi6s,2,quantile,c(.025,.975))
  qphi7<-apply(Phi7s,2,quantile,c(.025,.975))
  qphi8<-apply(Phi8s,2,quantile,c(.025,.975))
  
  #------------------------------------------------------------------------------#
  # Created County-Specific Trends #
  #--------------------------------#
  
  ts<--4:7
  num<-length(ts)
  
  MUPRED_NHW<-array(dim=c(lastit,length(ts),n))
  MUPRED_NHB<-array(dim=c(lastit,length(ts),n))
  
  DIFF_NHBW<-array(dim=c(lastit,length(ts),n))
  
  for (l in 1:lastit){
    alpha<-Alpha[l,]
    beta<-Beta[l,]
    phi1<-Phi1s[l,]
    phi2<-Phi2s[l,]
    phi3<-Phi3s[l,]
    phi4<-Phi4s[l,]
    phi5<-Phi5s[l,]
    phi6<-Phi6s[l,]
    phi7<-Phi7s[l,]
    phi8<-Phi8s[l,]
    for (k in 1:n){
      XNHW<-cbind(1,ts,0,0,rep(svis[k],num),0,0,0,0,0)
      XNHB<-cbind(1,ts,1,ts,rep(svis[k],num),0,0,0,0,0)
      
      # NHW
      eta1<-XNHW%*%alpha+phi1[k]+phi2[k]*ts
      mu1<-exp(eta1)/(1+exp(eta1))
      eta2<-XNHW%*%beta+phi5[k]+phi6[k]*ts
      mu2<-exp(eta2)/(1+exp(eta2))
      mu_nhw<-mu1*mu2*5
      
      MUPRED_NHW[l,,k]<-mu_nhw # Marginal Mean
      
      # NHB
      eta1<-XNHB%*%alpha+phi1[k]+phi2[k]*ts+phi3[k]+phi4[k]*ts
      mu1<-exp(eta1)/(1+exp(eta1))
      eta2<-XNHB%*%beta+phi5[k]+phi6[k]*ts+phi7[k]+phi8[k]*ts
      mu2<-exp(eta2)/(1+exp(eta2))
      mu_nhb<-mu1*mu2*5
      
      MUPRED_NHB[l,,k]<-mu_nhb # Marginal Mean
      
      # Difference in mean
      DIFF_NHBW[l,,k]<-mu_nhb-mu_nhw
    }
    if (l%%50==0) print(l)
  }
  
  # Posterior marginal CRFs for each county 
  mupred_nhw<-apply(MUPRED_NHW,3,colMeans)  # NHW
  mupred_nhb<-apply(MUPRED_NHB,3,colMeans)  # NHB
  mupred<-rbind(mupred_nhw,mupred_nhb)      # Combine
  # Posterior difference marginal CRFs
  diffpred_nhbw<-apply(DIFF_NHBW,3,colMeans)
  
  # 95 CrI - NHW
  mu_nhw_lower<-matrix(NA,num,n)
  for (l in 1:n) mu_nhw_lower[,l]<-apply(MUPRED_NHW[,,l],2,quantile,.025)
  mu_nhw_upper<-matrix(NA,num,n)
  for (l in 1:n) mu_nhw_upper[,l]<-apply(MUPRED_NHW[,,l],2,quantile,.975)
  # 95 CrI - NHB
  mu_nhb_lower<-matrix(NA,num,n)
  for (l in 1:n) mu_nhb_lower[,l]<-apply(MUPRED_NHB[,,l],2,quantile,.025)
  mu_nhb_upper<-matrix(NA,num,n)
  for (l in 1:n) mu_nhb_upper[,l]<-apply(MUPRED_NHB[,,l],2,quantile,.975)
  
  # 95 CrI - Difference
  diff_nhbw_lower<-matrix(NA,num,n)
  for (l in 1:n) diff_nhbw_lower[,l]<-apply(DIFF_NHBW[,,l],2,quantile,.025)
  diff_nhbw_upper<-matrix(NA,num,n)
  for (l in 1:n) diff_nhbw_upper[,l]<-apply(DIFF_NHBW[,,l],2,quantile,.975)
  
  mu_lower<-rbind(mu_nhw_lower,mu_nhb_lower)
  mu_upper<-rbind(mu_nhw_upper,mu_nhb_upper)
  
  
  dcounty<-data.frame(mu=c(mupred),
                      t=rep(rep(ts+5,2),n),
                      lb=c(mu_lower),
                      ub=c(mu_upper),
                      gp=rep(c(rep(1,num),rep(2,num)),n),
                      countyid=rep(1:n,each=num*2))
  
  dcounty$gp<-factor(dcounty$gp,levels=c(1,2),labels=c("NHW - Posterior Trend","NHB - Posterior Trend"))
  
  
  countylabel<-c("Abbeville","Aiken","Allendale","Anderson","Bamberg",
                 "Barnwell","Beaufort","Berkeley","Calhoun","Charleston",
                 "Cherokee","Chester","Chesterfield","Clarendon","Colleton",
                 "Darlington","Dillon","Dorchester","Edgefield","Fairfield",
                 "Florence","Georgetown","Greenville","Greenwood","Hampton",
                 "Horry","Jasper","Kershaw","Lancaster","Laurens",
                 "Lee","Lexington","McCormick","Marion","Marlboro",
                 "Newberry","Oconee","Orangeburg","Pickens","Richland",
                 "Saluda","Spartanburg","Sumter","Union","Williamsburg",
                 "York")
  
  dcounty$county<-factor(dcounty$countyid,levels=1:46,labels=countylabel)
  
  
  qrlab<-c("Q1, 2020", "Q2, 2020", "Q3, 2020", "Q4, 2020",
           "Q1, 2021", "Q2, 2021", "Q3, 2021", "Q4, 2021",
           "Q1, 2022", "Q2, 2022", "Q3, 2022", "Q4, 2022")
  
  # Selected counties
  plot1a=
    ggplot(dcounty[dcounty$county=="Chesterfield",],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Intervals",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Intervals"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,2),breaks = seq(0,2,by=0.5))+
    scale_color_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c( "solid", "solid","solid"),
                                  linewidth = c(1,1,3.5),
                                  shape=c(17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.75),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  plot1b=
    ggplot(dcounty[dcounty$county=="Clarendon",],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,1.25),breaks = seq(0,1.25,by=0.5))+
    scale_color_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c( "solid", "solid","solid"),
                                  linewidth = c(1,1,3.5),
                                  shape=c(17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.23),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  plot1c=
    ggplot(dcounty[dcounty$county=="Kershaw",],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,1.25),breaks = seq(0,1.25,by=0.5))+
    scale_color_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c( "solid", "solid","solid"),
                                  linewidth = c(1,1,3.5),
                                  shape=c(17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.23),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  plot1d=
    ggplot(dcounty[dcounty$county=="Jasper",],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,2),breaks = seq(0,2,by=0.5))+
    scale_color_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c( "solid", "solid","solid"),
                                  linewidth = c(1,1,3.5),
                                  shape=c(17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.75),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  
  ggarrange(plot1a+ rremove("xlab"),
            plot1b+ rremove("xlab")+rremove("ylab"),
            plot1c,
            plot1d+ rremove("ylab"),
            labels = c("(a)", "(b)", "(c)", "(d)"),
            font.label = list(size = 30),
            ncol = 2, nrow = 2,hjust = 0.01)
  
  # ggsave("",plot=pcounty,height = 16, width = 22,dpi = 300)
  
  #------------------------------------------------------------------------------#
  # Difference in mean CRF trend #
  #------------------------------#
  
  diffcounty<-data.frame(diff=c(diffpred_nhbw),
                         t=rep(ts+5,n),
                         lb=c(diff_nhbw_lower),
                         ub=c(diff_nhbw_upper),
                         countyid=rep(1:n,each=num))
  
  countylabel<-c("Abbeville","Aiken","Allendale","Anderson","Bamberg",
                 "Barnwell","Beaufort","Berkeley","Calhoun","Charleston",
                 "Cherokee","Chester","Chesterfield","Clarendon","Colleton",
                 "Darlington","Dillon","Dorchester","Edgefield","Fairfield",
                 "Florence","Georgetown","Greenville","Greenwood","Hampton",
                 "Horry","Jasper","Kershaw","Lancaster","Laurens",
                 "Lee","Lexington","McCormick","Marion","Marlboro",
                 "Newberry","Oconee","Orangeburg","Pickens","Richland",
                 "Saluda","Spartanburg","Sumter","Union","Williamsburg",
                 "York")
  
  diffcounty$county<-factor(diffcounty$countyid,levels=1:46,labels=countylabel)
  
  qrlab<-c("Q1, 2020", "Q2, 2020", "Q3, 2020", "Q4, 2020",
           "Q1, 2021", "Q2, 2021", "Q3, 2021", "Q4, 2021",
           "Q1, 2022", "Q2, 2022", "Q3, 2022", "Q4, 2022")
  
  # Selected counties 
  dplot1a=
    ggplot(diffcounty[diffcounty$county=="Chesterfield",],aes(x=t,y=diff))+
    geom_point(size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Intervals"),linetype="twodash",alpha=0.4,show.legend = T)+
    xlab("Yearly Quarter")+ylab("Difference in Mean CRF Score")+
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth=2)+
    scale_fill_manual(breaks = c("95% Credible Intervals"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(-0.2,1),breaks = seq(-0.2,1,by=0.2))+
    scale_color_manual(breaks = c("95% Credible Intervals"), 
                       values = c("grey45"))+
    guides(fill = guide_legend(title="Difference in mean CRF score",
                               override.aes = list(
                                 linetype = c("solid"),
                                 linewidth = c(3.5),
                                 shape=c(NA)),
                               reverse = F), color="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.85),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  
  dplot1b=
    ggplot(diffcounty[diffcounty$county=="Clarendon",],aes(x=t,y=diff))+
    geom_point(size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Intervals"),linetype="twodash",alpha=0.4,show.legend = T)+
    xlab("Yearly Quarter")+ylab("Difference in Mean CRF Score")+
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth=2)+
    scale_fill_manual(breaks = c("95% Credible Intervals"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(-0.2,1),breaks = seq(-0.2,1,by=0.2))+
    scale_color_manual(breaks = c("95% Credible Intervals"), 
                       values = c("grey45"))+
    guides(fill = guide_legend(title="Difference in mean CRF score",
                               override.aes = list(
                                 linetype = c("solid"),
                                 linewidth = c(3.5),
                                 shape=c(NA)),
                               reverse = F), color="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.85),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  dplot1c=
    ggplot(diffcounty[diffcounty$county=="Kershaw",],aes(x=t,y=diff))+
    geom_point(size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Intervals"),linetype="twodash",alpha=0.4,show.legend = T)+
    xlab("Yearly Quarter")+ylab("Difference in Mean CRF Score")+
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth=2)+
    scale_fill_manual(breaks = c("95% Credible Intervals"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(-0.2,1),breaks = seq(-0.2,1,by=0.2))+
    scale_color_manual(breaks = c("95% Credible Intervals"), 
                       values = c("grey45"))+
    guides(fill = guide_legend(title="Difference in mean CRF score",
                               override.aes = list(
                                 linetype = c("solid"),
                                 linewidth = c(3.5),
                                 shape=c(NA)),
                               reverse = F), color="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.85),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  dplot1d=
    ggplot(diffcounty[diffcounty$county=="Jasper",],aes(x=t,y=diff))+
    geom_point(size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Intervals"),linetype="twodash",alpha=0.4,show.legend = T)+
    xlab("Yearly Quarter")+ylab("Difference in Mean CRF Score")+
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth=2)+
    scale_fill_manual(breaks = c("95% Credible Intervals"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(-0.2,1),breaks = seq(-0.2,1,by=0.2))+
    scale_color_manual(breaks = c("95% Credible Intervals"), 
                       values = c("grey45"))+
    guides(fill = guide_legend(title="Difference in mean CRF score",
                               override.aes = list(
                                 linetype = c("solid"),
                                 linewidth = c(3.5),
                                 shape=c(NA)),
                               reverse = F), color="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.85),legend.text=element_text(size=20),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 35))+
    facet_wrap(~county)
  
  ggarrange(dplot1a+ rremove("xlab"),
            dplot1b+ rremove("xlab")+rremove("ylab"),
            dplot1c,
            dplot1d+ rremove("ylab"),
            labels = c("(a)", "(b)", "(c)", "(d)"),
            font.label = list(size = 30),
            ncol = 2, nrow = 2,hjust = 0.01)
  
  # ggsave(".jpg",height = 16, width = 22,dpi = 700)
  
  
  #------------------------------------------------------------------------------#
  # Predicted Mean County Over Time #
  #---------------------------------#
  
  library(berryFunctions)
  library(tmap)
  library(tigris)
  library(RColorBrewer)
  
  YPRED<-array(dim=c(lastit,N))
  for (l in 1:lastit){
    alpha<-Alpha[l,]
    beta<-Beta[l,]
    Phi1<-rep(Phi1s[l,],nis)
    Phi2<-rep(Phi2s[l,],nis)
    Phi3<-rep(Phi3s[l,],nis)
    Phi4<-rep(Phi4s[l,],nis)
    Phi5<-rep(Phi5s[l,],nis)
    Phi6<-rep(Phi6s[l,],nis)
    Phi7<-rep(Phi7s[l,],nis)
    Phi8<-rep(Phi8s[l,],nis)
    
    eta1<-X%*%alpha+Phi1+Phi2*t+Phi3*nhb+Phi4*int
    mu1<-exp(eta1)/(1+exp(eta1))
    eta2<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu2<-exp(eta2)/(1+exp(eta2))
    YPRED[l,]<-mu1*mu2*5
    
    if (l%%50==0) print(l)
  }
  
  # Average acorss whole maternal covariates
  ypred<-colMeans(YPRED)
  
  dat2<-data.frame(countyid=dsyn$countyid,
                  timeqr=dsyn$timeqr,
                  nhb=dsyn$nhb,
                  y=dsyn$y,
                  ypred=ypred)%>%
    filter(timeqr %in% c(1,4,8,12))%>%  # Only choose Q1,2020/ Q4,2020/ Q4,2021/ Q4,2022
    group_by(timeqr,nhb,countyid)%>%    # Should have 368 obs (4*2*46)
    summarise(mean_y=mean(y),             
              mean_ypred=mean(ypred))  # only 367 obs (No pregnancy record for NHW in county 33 at time 8, Q42021)
  
  
  # Making maps 
  # 46 (counties) * 4 (quarter times) * 2 (races)
  scfips<-c("45001","45003","45005","45007","45009",
            "45011","45013","45015","45017","45019","45021",
            "45023","45025","45027","45029","45031","45033",
            "45035","45037","45039","45041","45043","45045",
            "45047","45049","45051","45053","45055","45057",
            "45059","45061","45063","45065","45067","45069",
            "45071","45073","45075","45077","45079","45081",
            "45083","45085","45087","45089","45091")
  
  dat2$fips<-rep(scfips,8)
  dat2$GEOID<-as.character(dat2$fips)
  sc_counties_map0 <- counties(state = c('South Carolina'))
  sc_counties_map0$contfip<-as.numeric(sc_counties_map0$COUNTYFP)         # Create a numeric FIPS
  sc_counties_map1 <- sc_counties_map0[order(sc_counties_map0$contfip),]  # Sort dataset by FIPS: 001, 003, 005,...
  sc_counties_map8<-rbind(sc_counties_map1,sc_counties_map1,sc_counties_map1,sc_counties_map1, # Stack same dataset 8 times
                          sc_counties_map1,sc_counties_map1,sc_counties_map1,sc_counties_map1)
  sc_app_data <- cbind(sc_counties_map8,dat2)
  sc_app_data$nhb2<-factor(sc_app_data$nhb,levels = c(0,1),labels=c("NHW","NHB"))
  sc_app_data$timeqr2<-factor(sc_app_data$timeqr,levels = c(1,4,8,12),labels = c("Q1, 2020", "Q4, 2020","Q4, 2021", "Q4, 2022"))
  pal <- brewer.pal(5,"PuBu")  # specify the palette colors
  
  # Predicted mean CRFs
  tm_shape(sc_app_data)+
    tm_fill(c("mean_ypred"), midpoint = c(NA), title = c("Mean CRF Score"), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 2.5,
              legend.title.size = 4.5,
              # legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)+
    tm_facets(by=c("nhb2","timeqr2"),ncol=4)+
    tm_layout(panel.label.size = 4.5)
  
  
  #------------------------------------------------------------------------------#
