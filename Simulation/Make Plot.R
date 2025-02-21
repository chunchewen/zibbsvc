  # Make Figures
  # Import: MCMC Samples from Simulation-ZIBBSVC.R
  # File: Simulation-ZIBBSVC.Rda 
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
  set.seed(0425)
  n<-46 				                   # Number of county in SC
  nis<-c(599,2327,167,5688,338,575,3636,6555,311,11553,1594,615,912,741,1122,
         1969,954,4593,215,493,4048,1242,14043,1827,519,7375,567,1788,1098,1768,
         427,7614,97,857,366,957,1646,2205,2920,11013,397,9503,3264,693,734,2716)   # Number of pregnancy in each county
  ntrial<-5                        # Total number of CRFs (m, in the paper)
  id<-rep(1:n,nis)                 # County ID
  N<-length(id)                    # Number of total patients  
  
  # Covariate
  t<-sample(1:12,N,T)-5            # Centered at Q1, 2021
  nhb<-rbinom(N,1,0.40)            # NHB:NHW=0.4:0.6
  int<-t*nhb
  X<-cbind(1,t,nhb,int)            # fixed effect design matrix 
  p<-ncol(X)                       # N x p
  
  # Fixed Effect Parms
  true.alpha<-alpha<-c(0.35,0.40,0.20,-0.15)      # Binary component
  true.beta<-beta<-c(-0.75,0.15,-0.30,0.20)       # BB component 
  
  # Spatial random effect
  dir<-"..\\SC_agj.txt"            # Import SC adjacency matrix 
  A<-as.matrix(read.table(dir))    # SC Adjacency Table
  mi<-apply(A,1,sum)	             # No. neighbors
  Q<-as.spam(diag(mi))-as.spam(A) + diag(.0001,n) # Q=M-A
  
  cov<-matrix(c(0.30, 0.05, 0.10, 0.10, 0.05, 0.05, 0.05, 0.05, 
                0.05, 0.20, 0.05, 0.10, 0.05, 0.05, 0.05, -0.05,
                0.10, 0.05, 0.15, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.10, 0.10, 0.05, 0.20, 0.05, 0.05, 0.05, -0.05,
                0.05, 0.05, 0.05, 0.05, 0.25, 0.10, -0.10, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.10, 0.20, 0.05 , 0.05,
                0.05, 0.05, 0.05, 0.05, -0.10, 0.05, 0.20, 0.05,
                0.05,-0.05, 0.05,-0.05, 0.05, 0.05, 0.05 , 0.15),8,8)
  
  cov2cor(cov)                                # COV to CORRELATION
  covphi<-solve(Q)%x%cov			                # covariance of phis (8n x 8n)
  phi<-spam::rmvnorm(1,sigma=covphi)		      # 1 x 8n matrix of spatial effects
  phitmp<-matrix(phi, ncol=8, byrow=T)        # n x 8 matrix of spatial effects
  
  true.phi1<-phi1<-phitmp[,1]-mean(phitmp[,1])   # n x 1 phi1 vector -- Centered
  true.phi2<-phi2<-phitmp[,2]-mean(phitmp[,2])   # n x 1 phi2 vector, etc.
  true.phi3<-phi3<-phitmp[,3]-mean(phitmp[,3])  
  true.phi4<-phi4<-phitmp[,4]-mean(phitmp[,4])  
  true.phi5<-phi5<-phitmp[,5]-mean(phitmp[,5])  
  true.phi6<-phi6<-phitmp[,6]-mean(phitmp[,6])  
  true.phi7<-phi7<-phitmp[,7]-mean(phitmp[,7])  
  true.phi8<-phi8<-phitmp[,8]-mean(phitmp[,8])  
  
  true.phi<-cbind(true.phi1,true.phi2,true.phi3,true.phi4,true.phi5,true.phi6,true.phi7,true.phi8)
  
  Phi1<-rep(phi1,nis)
  Phi2<-rep(phi2,nis)
  Phi3<-rep(phi3,nis)
  Phi4<-rep(phi4,nis)
  Phi5<-rep(phi5,nis)
  Phi6<-rep(phi6,nis)
  Phi7<-rep(phi7,nis)
  Phi8<-rep(phi8,nis)
  
  # Overdispersion (0<rho<1)
  rho<-0.35     
  
  # Response
  # Binary component
  eta1<-X%*%alpha+Phi1+Phi2*t+Phi3*nhb+Phi4*int
  mu1<-exp(eta1)/(1+exp(eta1))
  u<-rbinom(N,1,mu1)                 # Abs. latent variable 
  N1<-sum(u)                         # true num of ppl at risk
  
  (N1/N)                             # proportion of ppl at risk
  (pstruct0<-1-mean(u))              # Proportion of structural zeros
  
  # BB component
  eta2<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
  mu2<-exp(eta2)/(1+exp(eta2))
  
  y<-rep(0,N)
  y[u==1]<-rbetabinom(N1,ntrial,mu2[u==1],rho)
  table(y)
  
  #---------------#
  # Store Samples #
  #---------------#
  nsim<-30000                   # Number of MCMC Iterations
  thin<-25				              # Thinnisng interval
  burn<-nsim/2   	              # Burnisn
  lastit<-(nsim-burn)/thin     	# Last stored value
  #------------------------------------------------------------------------------#
  
  #-------------------#
  # Load MCMC Samples #
  #-------------------#
  
  # dir.sav<-"...\\"
  # load(file=paste(dir.sav,"Simulation-ZIBBSVC.Rda",sep=""))           # Load MCMC samples
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
  
  #------------------------------------------------------------------------------#
  # Create Population-level Trends #
  #--------------------------------#
  num<-20                        # time grid 
  ts<-seq(-4,7,length.out=num)   # time variable
  
  XNHW<-cbind(1,ts,0,0)          # hypothetical NHW subject
  XNHB<-cbind(1,ts,1,ts)         # ... NHB subject
  
  # true trends b/w NHW and NHB 
  eta1_nhw<-XNHW%*%true.alpha;mu1_nhw<-exp(eta1_nhw)/(1+exp(eta1_nhw))
  eta2_nhw<-XNHW%*%true.beta;mu2_nhw<-exp(eta2_nhw)/(1+exp(eta2_nhw))
  true_mu_nhw<-mu1_nhw*mu2_nhw*5
  
  eta1_nhb<-XNHB%*%true.alpha;mu1_nhb<-exp(eta1_nhb)/(1+exp(eta1_nhb))
  eta2_nhb<-XNHB%*%true.beta;mu2_nhb<-exp(eta2_nhb)/(1+exp(eta2_nhb))
  true_mu_nhb<-mu1_nhb*mu2_nhb*5
  
  POP_MUPRED_NHW<-array(dim=c(lastit,length(ts)))
  POP_MUPRED_NHB<-array(dim=c(lastit,length(ts)))
  
  for (l in 1:lastit){
    alpha<-Alpha[l,]
    beta<-Beta[l,]
    
    # NHW
    eta1<-XNHW%*%alpha
    mu1<-exp(eta1)/(1+exp(eta1))    # E(R=1)
    eta2<-XNHW%*%beta
    mu2<-exp(eta2)/(1+exp(eta2))    # E(Y|R=1)
    mu_nhw<-mu1*mu2*5               # marginal mean E(Y)
    
    POP_MUPRED_NHW[l,]<-mu_nhw
    
    # NHB
    eta1<-XNHB%*%alpha
    mu1<-exp(eta1)/(1+exp(eta1))
    eta2<-XNHB%*%beta
    mu2<-exp(eta2)/(1+exp(eta2))
    mu_nhb<-mu1*mu2*5
    
    POP_MUPRED_NHB[l,]<-mu_nhb
    
    if (l%%50==0) print(l)
  }
  
  pop_mupred_nhw<-colMeans(POP_MUPRED_NHW)    # Population-level (average across counties)
  pop_mupred_nhb<-colMeans(POP_MUPRED_NHB)
  
  cri_pop_mupred_nhw<-apply(POP_MUPRED_NHW,2,quantile,c(.025,.975)) # CrIs
  cri_pop_mupred_nhb<-apply(POP_MUPRED_NHB,2,quantile,c(.025,.975))
  
  lower_pop_mupred_nhw<-cri_pop_mupred_nhw[1,] # LB: NHW
  upper_pop_mupred_nhw<-cri_pop_mupred_nhw[2,] # UB: NHW
  
  lower_pop_mupred_nhb<-cri_pop_mupred_nhb[1,] # LB: NHB
  upper_pop_mupred_nhb<-cri_pop_mupred_nhb[2,] # UB: NHB
  
  dplot_pop<-data.frame(mu=c(true_mu_nhw,true_mu_nhb,pop_mupred_nhw,pop_mupred_nhb),
                        lb=c(rep(NA,2*num),lower_pop_mupred_nhw,lower_pop_mupred_nhb),
                        ub=c(rep(NA,2*num),upper_pop_mupred_nhw,upper_pop_mupred_nhb),
                        t=rep(ts+5,4),
                        gp=rep(c(1,2,3,4),each=num))
  
  dplot_pop$gp<-factor(dplot_pop$gp,levels=c(1,2,3,4),labels=c("NHW - True Trend","NHB - True Trend",
                                                               "NHW - Posterior Trend","NHB - Posterior Trend"))
  
  qrlab<-c("Q1, 2020", "Q2, 2020", "Q3, 2020", "Q4, 2020",
           "Q1, 2021", "Q2, 2021", "Q3, 2021", "Q4, 2021",
           "Q1, 2022", "Q2, 2022", "Q3, 2022", "Q4, 2022")
  
  ggplot(dplot_pop,aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,label=qrlab)+scale_y_continuous(limits = c(0,5),breaks = seq(0,5,by=1))+
    scale_color_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("rosybrown2","royalblue","red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(15,16,17,18,NA))+
    guides(color = guide_legend(title="Population-level trends",
                                override.aes = list(
                                  linetype = c("solid", "solid", "solid", "solid","solid"),
                                  linewidth = c(1,1,1,1,3.5),
                                  shape=c(15,16,17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.20,0.75),legend.text=element_text(size=18),
          legend.title=element_text(size=20),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=18),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"))
  
  #---------------------------------------------------------------------#
  # Generate Map to compared simulated and estimated outcomes           #
  # Q1, 2020/ Q4, 2020/ Q4, 2021/ Q4, 2022                              #
  #---------------------------------------------------------------------#
  
  YPRED<-matrix(NA,lastit,N)
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
    
    # Estimated CRFs
    eta1<-X%*%alpha+rep(phi1,nis)+rep(phi2,nis)*t+rep(phi3,nis)*nhb+rep(phi4,nis)*int
    mu1<-exp(eta1)/(1+exp(eta1))
    eta2<-X%*%beta+rep(phi5,nis)+rep(phi6,nis)*t+rep(phi7,nis)*nhb+rep(phi8,nis)*int
    mu2<-exp(eta2)/(1+exp(eta2))
    nu<-mu1*mu2*5                             # E(Y)
    
    YPRED[l,]<-nu
    
    if (l%%50==0) print(l)
  }
  
  ypred<-colMeans(YPRED)                
  
  datmap<-data.frame(y=y,
                     ypred=ypred,
                     timeqr=t+5,
                     countyid=id)%>%
    filter(timeqr %in% c(1,4,8,12))%>%          # Only include Q1,2020/ Q4,2020/Q4,2021/Q4,2022
    group_by(countyid,timeqr)%>%                # group by county and quater time
    summarise(y=mean(y),ypred=mean(ypred))%>%   # averaage across county/quarter (true/estimated)
    pivot_wider(id_cols=countyid,
                names_from="timeqr",
                values_from=c("y","ypred"),
                names_glue="{timeqr}.{.value}")%>%
    pivot_longer(cols=`1.y`:`12.ypred`,
                 names_to="Time")
  # county SC FIPS
  scfips<-c("45001","45003","45005","45007","45009",
            "45011","45013","45015","45017","45019","45021",
            "45023","45025","45027","45029","45031","45033",
            "45035","45037","45039","45041","45043","45045",
            "45047","45049","45051","45053","45055","45057",
            "45059","45061","45063","45065","45067","45069",
            "45071","45073","45075","45077","45079","45081",
            "45083","45085","45087","45089","45091")
  
  
  datmap$fips<-rep(scfips,each=8)                                          # 4 quarter times for true map and 4 for estiamted map
  datmap$GEOID<-as.character(datmap$fips)
  sc_counties_map0 <- counties(state = c('South Carolina'))                # Extract SC info
  sc_counties_map0$contfip<-as.numeric(sc_counties_map0$COUNTYFP)          # Create a numeric FIPS
  sc_counties_map1 <- sc_counties_map0[order(sc_counties_map0$contfip),]   # Sort dataset by FIPS: 001, 003, 005,...
  sc_counties_map8<-rbind(sc_counties_map1,sc_counties_map1,sc_counties_map1,sc_counties_map1, # Stack same dataset 8 times
                          sc_counties_map1,sc_counties_map1,sc_counties_map1,sc_counties_map1)
  sc_counties_map8_v2<-sc_counties_map8[order(sc_counties_map8$contfip),]  # Sort FIPS again 001x8, 003x8, 005x8
  sc_counties_map8_v2$group<-rep(c("Q1, 2020","Q4, 2020","Q4, 2021", "Q4, 2022"),46*2)
  sc_counties_map8_v2$datgen<-rep(c(rep("Simulated",4),rep("Estimated",4)),46)
  sc_counties_map8_v2$datgen2<-factor(sc_counties_map8_v2$datgen,levels=c("Simulated","Estimated"))
  sc_sim_data <- cbind(sc_counties_map8_v2,datmap)
  
  pal <- brewer.pal(5,"PuBu")  # specify the palette colors
  
  pmap=tm_shape(sc_sim_data)+
    tm_fill(c("value"), midpoint = c(NA), title = c("Mean CRF Score"), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 2.5,
              legend.title.size = 4.5,
              main.title.fontface = "bold",
              main.title.position = "center",
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)+
    tm_facets(by=c("datgen2","group"),ncol=4)+
    tm_layout(panel.label.size = 4.5)
  # tmap_save(pmap, filename = "xxx.png", width = 16, height = 10, dpi = 300)  # save high resolution map
  
  #------------------------------------------------------------------------------#
  # Create Map Spatially Random Effects       #
  # compare true vs estiamted spatial effects #
  #-------------------------------------------#
  mphi1<-colMeans(Phi1s)
  mphi2<-colMeans(Phi2s)
  mphi3<-colMeans(Phi3s)
  mphi4<-colMeans(Phi4s)
  mphi5<-colMeans(Phi5s)
  mphi6<-colMeans(Phi6s)
  mphi7<-colMeans(Phi7s)
  mphi8<-colMeans(Phi8s)
  
  dat<-data.frame(y=y,
                  countyid=id)
  
  dmap<-dat %>% group_by(countyid)%>%
    summarize(n=n(),
              mMRF=mean(y))
  
  
  # true random effects
  dmap$phi1<-true.phi1
  dmap$phi2<-true.phi2
  dmap$phi3<-true.phi3
  dmap$phi4<-true.phi4
  dmap$phi5<-true.phi5
  dmap$phi6<-true.phi6
  dmap$phi7<-true.phi7
  dmap$phi8<-true.phi8
  
  dmap$phi1p<-mphi1
  dmap$phi2p<-mphi2
  dmap$phi3p<-mphi3
  dmap$phi4p<-mphi4
  dmap$phi5p<-mphi5
  dmap$phi6p<-mphi6
  dmap$phi7p<-mphi7
  dmap$phi8p<-mphi8
  
  scfips<-c("45001","45003","45005","45007","45009",
            "45011","45013","45015","45017","45019","45021",
            "45023","45025","45027","45029","45031","45033",
            "45035","45037","45039","45041","45043","45045",
            "45047","45049","45051","45053","45055","45057",
            "45059","45061","45063","45065","45067","45069",
            "45071","45073","45075","45077","45079","45081",
            "45083","45085","45087","45089","45091")
  
  
  dmap$fips<-scfips
  dmap$GEOID<-as.character(dmap$fips)
  sc_counties_map <- counties(state = c('South Carolina'))
  sc_mcmf_data <- left_join(sc_counties_map, dmap, by = 'GEOID')
  pal <- brewer.pal(5,"PuBu")                  # specify the palette colors

  
  # Combine random effects for phi4 and phi8 (for paper) 
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi4","phi4p","phi8","phi8p"), midpoint = c(NA), title =  c(""), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              panel.labels = c(expression(paste("(a) Simulated ", phi[1*i*4], sep = " ")), 
                               expression(paste("(b) Estimated ", phi[1*i*4], sep = " ")),
                               expression(paste("(c) Simulated ", phi[2*i*4], sep = " ")), 
                               expression(paste("(d) Estimated ", phi[2*i*4], sep = " "))),
              panel.label.bg.color = c('white'),
              panel.label.fontface = c('bold'),
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.position = c(0.01, 0.10),
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Combine random effects for phi1 to phi8 (for supplement)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi1","phi1p","phi2","phi2p","phi3","phi3p","phi4","phi4p",
              "phi5","phi5p","phi6","phi6p","phi7","phi7p","phi8","phi8p"), midpoint = c(NA), title =  c(""), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              panel.labels = c(expression(paste("(A) Simulated ", phi[1*i*1], sep = " ")), 
                               expression(paste("(a) Estimated ", phi[1*i*1], sep = " ")),
                               expression(paste("(B) Simulated ", phi[1*i*2], sep = " ")), 
                               expression(paste("(b) Estimated ", phi[1*i*2], sep = " ")),
                               expression(paste("(C) Simulated ", phi[1*i*3], sep = " ")), 
                               expression(paste("(c) Estimated ", phi[1*i*3], sep = " ")),
                               expression(paste("(D) Simulated ", phi[1*i*4], sep = " ")), 
                               expression(paste("(d) Estimated ", phi[1*i*4], sep = " ")),
                               expression(paste("(E) Simulated ", phi[2*i*1], sep = " ")), 
                               expression(paste("(e) Estimated ", phi[2*i*1], sep = " ")),
                               expression(paste("(F) Simulated ", phi[2*i*2], sep = " ")), 
                               expression(paste("(f) Estimated ", phi[2*i*2], sep = " ")),
                               expression(paste("(G) Simulated ", phi[2*i*3], sep = " ")), 
                               expression(paste("(g) Estimated ", phi[2*i*3], sep = " ")),
                               expression(paste("(H) Simulated ", phi[2*i*4], sep = " ")), 
                               expression(paste("(h) Estimated ", phi[2*i*4], sep = " "))),
              panel.label.bg.color = c('white'),
              panel.label.fontface = c('bold'),
              legend.outside.size=0.25,
              legend.text.size = .75,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.position = c(0.01, 0.10),
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # tmap_save(prand, filename = "xxx.png", width = 16, height = 10, dpi = 300)  # save high resolution map
  
  
  #------------------------------------------------------------------------------#
  # County Specific Trends #
  #------------------------#
  
  num<-12                           # time grid
  ts<-seq(-4,7,length.out=num)      # time period
  
  XNHW_county<-kronecker(rep(1,n),cbind(1,ts,0,0))  # NHW design matrix for each county 
  XNHB_county<-kronecker(rep(1,n),cbind(1,ts,1,ts)) # NHB design matrix for each county 
  
  num_s<-rep(num,n)
  
  # county-specific true trends: NHW 
  eta1_nhw_county<-XNHW_county%*%true.alpha+rep(true.phi1,num_s)+rep(true.phi2,num_s)*ts
  mu1_nhw_county<-exp(eta1_nhw_county)/(1+exp(eta1_nhw_county))
  eta2_nhw_county<-XNHW_county%*%true.beta+rep(true.phi5,num_s)+rep(true.phi6,num_s)*ts
  mu2_nhw_county<-exp(eta2_nhw_county)/(1+exp(eta2_nhw_county))
  true_mu_nhw_county<-matrix(mu1_nhw_county*mu2_nhw_county*5,nrow=num,ncol=n)
  # county-specific true trends: NHB
  eta1_nhb_county<-XNHB_county%*%true.alpha+rep(true.phi1,num_s)+rep(true.phi2,num_s)*ts+rep(true.phi3,num_s)+rep(true.phi4,num_s)*ts
  mu1_nhb_county<-exp(eta1_nhb_county)/(1+exp(eta1_nhb_county))
  eta2_nhb_county<-XNHB_county%*%true.beta+rep(true.phi5,num_s)+rep(true.phi6,num_s)*ts+rep(true.phi7,num_s)+rep(true.phi8,num_s)*ts
  mu2_nhb_county<-exp(eta2_nhb_county)/(1+exp(eta2_nhb_county))
  true_mu_nhb_county<-matrix(mu1_nhb_county*mu2_nhb_county*5,nrow=num,ncol=n)
  
  
  # county-specific estiamted trends: NHW vs NHB
  true_mu<-rbind(true_mu_nhw_county,true_mu_nhb_county)
  
  MUPRED_NHW<-array(dim=c(lastit,length(ts),n))
  MUPRED_NHB<-array(dim=c(lastit,length(ts),n))
  
  
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
      XNHW<-cbind(1,ts,0,0)
      XNHB<-cbind(1,ts,1,ts)
      
      # NHW
      eta1<-XNHW%*%alpha+phi1[k]+phi2[k]*ts
      mu1<-exp(eta1)/(1+exp(eta1))
      eta2<-XNHW%*%beta+phi5[k]+phi6[k]*ts
      mu2<-exp(eta2)/(1+exp(eta2))
      mu_nhw<-mu1*mu2*5
      
      MUPRED_NHW[l,,k]<-mu_nhw
      
      # NHB
      eta1<-XNHB%*%alpha+phi1[k]+phi2[k]*ts+phi3[k]+phi4[k]*ts
      mu1<-exp(eta1)/(1+exp(eta1))
      eta2<-XNHB%*%beta+phi5[k]+phi6[k]*ts+phi7[k]+phi8[k]*ts
      mu2<-exp(eta2)/(1+exp(eta2))
      mu_nhb<-mu1*mu2*5
      
      MUPRED_NHB[l,,k]<-mu_nhb
    }
    if (l%%50==0) print(l)
  }
  
  
  mupred_nhw<-apply(MUPRED_NHW,3,colMeans)
  mupred_nhb<-apply(MUPRED_NHB,3,colMeans)
  mupred<-rbind(mupred_nhw,mupred_nhb)
  
  
  apply(MUPRED_NHW[,,1],2,quantile,c(.025,.975)) # countyid 
  
  mu_nhw_lower<-matrix(NA,num,n)
  for (l in 1:n) mu_nhw_lower[,l]<-apply(MUPRED_NHW[,,l],2,quantile,.025)
  mu_nhw_upper<-matrix(NA,num,n)
  for (l in 1:n) mu_nhw_upper[,l]<-apply(MUPRED_NHW[,,l],2,quantile,.975)
  
  mu_nhb_lower<-matrix(NA,num,n)
  for (l in 1:n) mu_nhb_lower[,l]<-apply(MUPRED_NHB[,,l],2,quantile,.025)
  mu_nhb_upper<-matrix(NA,num,n)
  for (l in 1:n) mu_nhb_upper[,l]<-apply(MUPRED_NHB[,,l],2,quantile,.975)
  
  mu_lower<-rbind(mu_nhw_lower,mu_nhb_lower)
  mu_upper<-rbind(mu_nhw_upper,mu_nhb_upper)
  
  
  dcounty<-data.frame(mu=c(c(true_mu),c(mupred)),
                      t=rep(ts+5,n*4),
                      lb=c(rep(NA,num*n*2),c(mu_lower)),
                      ub=c(rep(NA,num*n*2),c(mu_upper)),
                      gp=c(rep(c(rep(1,num),rep(2,num)),n),rep(c(rep(3,num),rep(4,num)),n)),
                      countyid=rep(rep(1:n,each=2*num),2))
  
  dcounty$gp<-factor(dcounty$gp,levels=c(1,2,3,4),labels=c("NHW - True Trend","NHB - True Trend",
                                                           "NHW - Posterior Trend","NHB - Posterior Trend"))
  

  
  qrlab<-c("Q1, 2020", "Q2, 2020", "Q3, 2020", "Q4, 2020",
           "Q1, 2021", "Q2, 2021", "Q3, 2021", "Q4, 2021",
           "Q1, 2022", "Q2, 2022", "Q3, 2022", "Q4, 2022")
  
  plot1a=
    ggplot(dcounty[dcounty$countyid==18,],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,5),breaks = seq(0,5,by=1))+
    scale_color_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("rosybrown2","royalblue","red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(15,16,17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c("solid", "solid", "solid", "solid","solid"),
                                  linewidth = c(1,1,1,1,3.5),
                                  shape=c(15,16,17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.75),legend.text=element_text(size=22),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=18),
          axis.title=element_text(size=22),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 30))+
    facet_wrap(~countyid)
  
  plot1b=
    ggplot(dcounty[dcounty$countyid==28,],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,5),breaks = seq(0,5,by=1))+
    scale_color_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("rosybrown2","royalblue","red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(15,16,17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c("solid", "solid", "solid", "solid","solid"),
                                  linewidth = c(1,1,1,1,3.5),
                                  shape=c(15,16,17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.75),legend.text=element_text(size=22),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=18),
          axis.title=element_text(size=22),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 30))+
    facet_wrap(~countyid)
  
  plot1c=
    ggplot(dcounty[dcounty$countyid==32,],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,5),breaks = seq(0,5,by=1))+
    scale_color_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("rosybrown2","royalblue","red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(15,16,17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c("solid", "solid", "solid", "solid","solid"),
                                  linewidth = c(1,1,1,1,3.5),
                                  shape=c(15,16,17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.70),legend.text=element_text(size=22),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=18),
          axis.title=element_text(size=22),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 30))+
    facet_wrap(~countyid)
  
  plot1d=
    ggplot(dcounty[dcounty$countyid==38,],aes(x=t,y=mu,col=gp))+
    geom_point(aes(shape=gp),size=4)+
    geom_line(aes(linetype="solid"),linewidth=1)+
    geom_ribbon(aes(ymin = lb, ymax = ub,col="95% Credible Intervals",fill="95% Credible Interval",group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    xlab("Yearly Quarter")+ylab("Mean CRF Score")+
    scale_fill_manual(breaks = c("95% Credible Interval"), 
                      values = c("grey36"))+
    scale_x_continuous(breaks = 1:12,labels=qrlab)+scale_y_continuous(limits = c(0,5),breaks = seq(0,5,by=1))+
    scale_color_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c("rosybrown2","royalblue","red4","blue4","grey45"))+
    scale_shape_manual(breaks = c("NHW - True Trend","NHB - True Trend",
                                  "NHW - Posterior Trend","NHB - Posterior Trend","95% Credible Intervals"), 
                       values = c(15,16,17,18,NA))+
    guides(color = guide_legend(title="Race-Ethnicity Group",
                                override.aes = list(
                                  linetype = c("solid", "solid", "solid", "solid","solid"),
                                  linewidth = c(1,1,1,1,3.5),
                                  shape=c(15,16,17,18,NA)),
                                reverse = F), fill="none",shape="none",linetype="none",group="none")+
    theme_gray(base_size = 14)+
    theme(legend.key = element_rect(fill = "white"),
          legend.key.width = unit(10,"mm"),
          legend.position = c(0.25,0.70),legend.text=element_text(size=22),
          legend.title=element_text(size=21),
          axis.text.x = element_text(angle = 70, hjust=1),
          axis.text=element_text(size=18),
          axis.title=element_text(size=22),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"),
          strip.text = element_text(size = 30))+
    facet_wrap(~countyid)
  
  ggarrange(plot1a+ rremove("xlab"),
            plot1b+ rremove("xlab")+rremove("ylab"),
            plot1c,
            plot1d+ rremove("ylab"),
            labels = c("(a)", "(b)", "(c)", "(d)"),
            font.label = list(size = 25),
            ncol = 2, nrow = 2,hjust = 0.0)
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  
