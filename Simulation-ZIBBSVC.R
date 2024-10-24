  # Zero-Inflated Beta-Binomial Spatially Varying Coefficients (ZIBB-SVC) Model
  # Spatially varying space x time in both component
  # Simulation Analysis 
  # Date: 10-22-2024
  #------------------------------------------------------------------------------#
  library(aod)
  library(VGAM)       # For rbetabinom()
  library(mcmcplots)
  library(spam)
  library(mvtnorm)    # For rmvt()
  library(msm)        # For rtnorm()
  library(BayesLogit) # For rpg() 
  library(MCMCpack)   # For riwsh()
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
  dir<-"C:\\Users\\chech\\OneDrive - Medical University of South Carolina\\Research\\ZIBBST Model\\R Code\\Personal\\Simulation\\SC_adj.txt"
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
  
  #-----------------------------#
  # Histogram of simulated CRFs #
  #-----------------------------#
  pbar<-data.frame(y=c(table(y)/sum(table(y))),
                   score=c(0,1,2,3,4,5))
  
  ggplot(pbar, aes(x=score, y=y)) + 
    geom_bar(stat = "identity",fill="#006633")+
    xlab("Simulated Cardiometabolic Risk Factor Score")+ylab("Percentage (%)")+
    scale_x_continuous(breaks = 0:5)+
    scale_y_continuous(breaks = seq(0,0.6,by=0.1))+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"))   
  #------------------------------------------------------------------------------#
  # MCMC prep #
  #-----------#
  # Priors
  alpha0<-rep(0,p)
  T0a<-diag(0.01,p)        # prior precision for alpha
  beta0<-rep(0,p)
  V0b<-diag(100,p)         # prior variance for beta
  c0<-9                    # prior for Wishart(nu0=diag(q),C0=q+1)
  d0<-diag(8)
  
  # Init
  alpha<-rep(0,p)
  beta<-rep(0,p)
  rho<-0.5
  A1<-A2<-0               # Acceptance counter (beta,rho)
  Accu1<-rep(0,n)         # Acceptance counter (phi5, each county)
  Accu2<-rep(0,n)         # Acceptance counter (phi6, each county)
  Accu3<-rep(0,n)         # Acceptance counter (phi7, each county)
  Accu4<-rep(0,n)         # Acceptance counter (phi8, each county)
  
  # proposal
  covb<-diag(.05,p)       # proposal var for beta
  sigmar0<-0.0001         # proposal var for rho
  s5<-0.65                # .. for phi5
  s6<-0.25                # .. for phi6...etc
  s7<-0.45
  s8<-0.30
  
  y1<-rbinom(N,1,.5)        # At risk indicator
  y1[y>0]<-1                # If y>0, then patient is at risk w.p. 1
  n0<-length(y[y==0])       # Number of observed 0's
  
  phi_init<-spam::rmvnorm(1, sigma=diag(.1, 8*n))	  # Random effects
  phi_init<-matrix(phi_init, ncol=8, byrow=T)       # n x 8 matrix of spatial effects
  phi1<-phi_init[,1]-mean(phi_init[,1])             # centering...
  phi2<-phi_init[,2]-mean(phi_init[,2])
  phi3<-phi_init[,3]-mean(phi_init[,3])
  phi4<-phi_init[,4]-mean(phi_init[,4])
  phi5<-phi_init[,5]-mean(phi_init[,5])
  phi6<-phi_init[,6]-mean(phi_init[,6])
  phi7<-phi_init[,7]-mean(phi_init[,7])
  phi8<-phi_init[,8]-mean(phi_init[,8])
  
  phimat<-cbind(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8)
  
  Phi1<-rep(phi1,nis)
  Phi2<-rep(phi2,nis)
  Phi3<-rep(phi3,nis)
  Phi4<-rep(phi4,nis)
  Phi5<-rep(phi5,nis)
  Phi6<-rep(phi6,nis)
  Phi7<-rep(phi7,nis)
  Phi8<-rep(phi8,nis)
  
  phibar<-matrix(0,n,8)                # Cond. of MCAR (mean of random effects for neighboring counties)
  for (j in 1:n){
    phibar[j,1]<-mean(phi1[which(A[j,]==1)])
    phibar[j,2]<-mean(phi2[which(A[j,]==1)])
    phibar[j,3]<-mean(phi3[which(A[j,]==1)])
    phibar[j,4]<-mean(phi4[which(A[j,]==1)])
    phibar[j,5]<-mean(phi5[which(A[j,]==1)])
    phibar[j,6]<-mean(phi6[which(A[j,]==1)])
    phibar[j,7]<-mean(phi7[which(A[j,]==1)])
    phibar[j,8]<-mean(phi8[which(A[j,]==1)])
  }
  
  # Random variance
  sigmaphi<-cov(phimat)
  
  #---------------#
  # Store Samples #
  #---------------#
  nsim<-30000                   # Number of MCMC Iterations
  thin<-25				              # Thinnisng interval
  burn<-nsim/2   	              # Burnisn
  lastit<-(nsim-burn)/thin     	# Last stored value
  Alpha<-matrix(NA,lastit,p)
  Betatmp<-matrix(NA,nsim,p)
  Beta<-matrix(NA,lastit,p)
  Rho<-rep(NA,lastit)
  S2Phi<-matrix(0,lastit,64)
  Phi1s<-matrix(NA,lastit,n)
  Phi2s<-matrix(NA,lastit,n)
  Phi3s<-matrix(NA,lastit,n)
  Phi4s<-matrix(NA,lastit,n)
  Phi5s<-matrix(NA,lastit,n)
  Phi6s<-matrix(NA,lastit,n)
  Phi7s<-matrix(NA,lastit,n)
  Phi8s<-matrix(NA,lastit,n)
  
  set.seed(1234)
  time.start<-proc.time()
  for (i in 1:nsim){
    #----------------#
    # Binary - Fixed #
    #----------------#
    mu<-X%*%alpha+Phi1+Phi2*t+Phi3*nhb+Phi4*int
    omega<-rpg(N,1,mu)
    z<-(y1-1/2)/omega
    v<-solve(T0a+crossprod(X*sqrt(omega)))
    m<-v%*%(T0a%*%alpha0+t(sqrt(omega)*X)%*%c(sqrt(omega)*(z-Phi1-Phi2*t-Phi3*nhb-Phi4*int)))
    alpha<-c(rmvnorm(1,m,v))
    
    #---------------#
    # Binary - phi1 #
    #---------------#
    priorprec<-c(1/(sigmaphi[1,1]-sigmaphi[1,-1]%*%solve(sigmaphi[-1,-1])%*%sigmaphi[-1,1]))*Q  # cond. precision phi1|rest 
    priormean<-phimat[,-1]%*%t(sigmaphi[1,-1]%*%solve(sigmaphi[-1,-1]))                         # cond. mean (n x 1)   
    prec<-priorprec+as.spam(diag(tapply(omega,id,sum),n,n))
    mb<-c(priorprec%*%priormean)+tapply(omega*(z-X%*%alpha-Phi2*t-Phi3*nhb-Phi4*int),id,sum)
    phi1<-rmvnorm.canonical(1, mb, prec)[1,]
    
    # Center phi1 and update phi1bar
    phi1<-phi1-mean(phi1)
    Phi1<-rep(phi1,nis)
    for (j in 1:n) phibar[j,1]<-mean(phi1[which(A[j,]==1)])
    phimat[,1]<-phi1
    
    #---------------#
    # Binary - phi2 #
    #---------------#
    priorprec<-c(1/(sigmaphi[2,2]-sigmaphi[2,-2]%*%solve(sigmaphi[-2,-2])%*%sigmaphi[-2,2]))*Q # cond. precision phi2|rest 
    priormean<-phimat[,-2]%*%t(sigmaphi[2,-2]%*%solve(sigmaphi[-2,-2]))                        # cond. mean (n x 1)
    prec<-priorprec+as.spam(diag(tapply(omega*t^2,id,sum),n,n))
    mb<-c(priorprec%*%priormean)+tapply(omega*t*(z-X%*%alpha-Phi1-Phi3*nhb-Phi4*int),id,sum)
    phi2<-rmvnorm.canonical(1, mb, prec)[1,]
    
    # Center phi2 and update phi2bar
    phi2<-phi2-mean(phi2)
    Phi2<-rep(phi2,nis)
    for (j in 1:n) phibar[j,2]<-mean(phi2[which(A[j,]==1)])
    phimat[,2]<-phi2
    
    #---------------#
    # Binary - phi3 #
    #---------------#
    priorprec<-c(1/(sigmaphi[3,3]-sigmaphi[3,-3]%*%solve(sigmaphi[-3,-3])%*%sigmaphi[-3,3]))*Q  # cond. precision phi3|rest 
    priormean<-phimat[,-3]%*%t(sigmaphi[3,-3]%*%solve(sigmaphi[-3,-3]))                         # cond. mean (n x 1)
    prec<-priorprec+as.spam(diag(tapply(omega*nhb^2,id,sum),n,n))
    mb<-c(priorprec%*%priormean)+tapply(omega*nhb*(z-X%*%alpha-Phi1-Phi2*t-Phi4*int),id,sum)
    phi3<-rmvnorm.canonical(1, mb, prec)[1,]
    
    # Center phi3 and update phi3bar
    phi3<-phi3-mean(phi3)
    Phi3<-rep(phi3,nis)
    for (j in 1:n) phibar[j,3]<-mean(phi3[which(A[j,]==1)])
    phimat[,3]<-phi3
    
    #---------------#
    # Binary - phi4 #
    #---------------#
    priorprec<-c(1/(sigmaphi[4,4]-sigmaphi[4,-4]%*%solve(sigmaphi[-4,-4])%*%sigmaphi[-4,4]))*Q # cond. precision ph41|rest 
    priormean<-phimat[,-4]%*%t(sigmaphi[4,-4]%*%solve(sigmaphi[-4,-4]))                        # cond. mean (n x 1)
    prec<-priorprec+as.spam(diag(tapply(omega*int^2,id,sum),n,n))
    mb<-c(priorprec%*%priormean)+tapply(omega*int*(z-X%*%alpha-Phi1-Phi2*t-Phi3*nhb),id,sum)
    phi4<-rmvnorm.canonical(1, mb, prec)[1,]
    
    # Center phi4 and update phi4bar
    phi4<-phi4-mean(phi4)
    Phi4<-rep(phi4,nis)
    for (j in 1:n) phibar[j,4]<-mean(phi4[which(A[j,]==1)])
    phimat[,4]<-phi4
    
    #--------------------#
    # Update at-risk ind #
    #--------------------#
    eta1<-X%*%alpha+Phi1+Phi2*t+Phi3*nhb+Phi4*int
    eta2<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    pi<-pmax(0.001,pmin(0.999,1/(1+exp(-eta1))))    # at-risk probability
    pr<-1/(1+exp(-eta2))
    q<-dbetabinom(0,size=ntrial,prob=pr,rho=rho)    
    theta<-pi*q/(pi*q+1-pi)                         # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
    y1[y==0]<-rbinom(n0,1,theta[y==0])              # If y=0, then draw a "chance zero" w.p. theta, if y=1, then y1=1
    nis1<-tapply(y1,id,sum)
    
    #-------------#
    #  BB - fixed #
    #-------------#
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lold<-sum(dbetabinom(y[y1==1],size=ntrial,prob=mu[y1==1],rho=rho,log=T))
    
    betanew<-beta+rmvnorm(1,sigma=.05*covb)        # Draw from "symmetric" MV dist
    eta<-X%*%c(betanew)+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lnew<-sum(dbetabinom(y[y1==1],size=ntrial,prob=mu[y1==1],rho=rho,log=T))
    
    
    # Acceptance prob on log scale =log(lnew x prior) - log (lold x prior)
    ratio<-lnew+dmvnorm(betanew,beta0,V0b,log=T)-(lold+dmvnorm(beta,beta0,V0b,log=T))
    if(log(runif(1))<ratio) {
      beta<-c(betanew)
      A1<-A1+1
    }
    
    Betatmp[i,]<-beta
    if (i==nsim/2) covb<-cov(Betatmp[(nsim/4+1):nsim/2,])  # Update proposal cov
    
    #--------#
    # BB rho #
    #--------#
    # Current likelihood
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lold<-sum(dbetabinom(y[y1==1],size=ntrial,prob=mu[y1==1],rho=rho,log=T))
    
    # Draw candidate rho and compute likelihood from truncated noraml
    rhonew<-rtnorm(1,rho,sqrt(sigmar0),0,1)           # Draw from truncated normal
    lnew<-sum(dbetabinom(y[y1==1],size=ntrial,prob=mu[y1==1],rho=rhonew,log=T))
    
    # Acceptance prob on log scale
    ratio<-lnew-lold+dtnorm(rho,rhonew,sqrt(sigmar0),0,1,log=T)-dtnorm(rhonew,rho,sqrt(sigmar0),0,1,log=T)
    if(log(runif(1))<ratio) {
      rho<-rhonew
      A2<-A2+1         
    }
    
    #---------#
    # BB phi5 #
    #---------#
    phi5new<-phi5+s5*rt(n, 3)         # Proposal distribution
    eta<-X%*%beta+rep(phi5new,nis)+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lnew<-rep(0,n)                    # account for empty counties - their likelihood contribution is 0
    lnew[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lold<-rep(0,n)                    # account for empty counties - their likelihood contribution is 0
    lold[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    for (j in 1:n){
      phibar[j,5]<-mean(phi5[which(A[j,]==1)])  # update phi5bar[j] before updating phi5[j]
      m<-phibar[j,5]+sigmaphi[5,-5]%*%solve(sigmaphi[-5,-5])%*%(phimat[j,-5]-phibar[j,-5]) # cond. mean phi5i|rest
      v<-(sigmaphi[5,5]-sigmaphi[5,-5]%*%solve(sigmaphi[-5,-5])%*%sigmaphi[-5,5])/mi[j]    # cond. var
      pnew<-dnorm(phi5new[j],m, sqrt(v), log=T)
      pold<-dnorm(phi5[j],m, sqrt(v), log=T)
      ratio<-(lnew[j]+pnew)-(lold[j]+pold)
      if (log(runif(1))<ratio) {
        phi5[j]<-phi5new[j]
        Accu1[j]<-Accu1[j]+1                    # keep track of acceptance of phi5 for each county
      }
    }
    
    # Center phi5
    phi5<-phi5-mean(phi5)
    Phi5<-rep(phi5,nis)
    phimat[,5]<-phi5
    
    #---------#
    # BB phi6 #
    #---------#
    phi6new<-phi6+s6*rt(n, 3)                   # Proposal distribution
    eta<-X%*%beta+Phi5+rep(phi6new,nis)*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lnew<-rep(0,n)                              # account for empty counties - their likelihood contribution is 0
    lnew[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lold<-rep(0,n)                              # account for empty counties - their likelihood contribution is 0
    lold[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    for (j in 1:n){
      phibar[j,6]<-mean(phi6[which(A[j,]==1)])  # update phi6bar[j] before updating phi6[j]
      m<-phibar[j,6]+sigmaphi[6,-6]%*%solve(sigmaphi[-6,-6])%*%(phimat[j,-6]-phibar[j,-6]) # cond. mean phi7i|rest
      v<-(sigmaphi[6,6]-sigmaphi[6,-6]%*%solve(sigmaphi[-6,-6])%*%sigmaphi[-6,6])/mi[j]    # cond. var
      pnew<-dnorm(phi6new[j],m, sqrt(v), log=T)
      pold<-dnorm(phi6[j],m, sqrt(v), log=T)
      ratio<-(lnew[j]+pnew)-(lold[j]+pold)
      if (log(runif(1))<ratio) {
        phi6[j]<-phi6new[j]
        Accu2[j]<-Accu2[j]+1                    # keep track of acceptance of phi for each county
      }
    }
    
    # Center phi6
    phi6<-phi6-mean(phi6)
    Phi6<-rep(phi6,nis)
    phimat[,6]<-phi6
    
    #---------#
    # BB phi7 #
    #---------#
    phi7new<-phi7+s7*rt(n, 3)        # Proposal distribution
    eta<-X%*%beta+Phi5+Phi6*t+rep(phi7new,nis)*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lnew<-rep(0,n)                 # account for empty counties - their likelihood contribution is 0
    lnew[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lold<-rep(0,n)                    # account for empty counties - their likelihood contribution is 0
    lold[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    for (j in 1:n){
      phibar[j,7]<-mean(phi7[which(A[j,]==1)])  # update phi7bar[j] before updating phi7[j]
      m<-phibar[j,7]+sigmaphi[7,-7]%*%solve(sigmaphi[-7,-7])%*%(phimat[j,-7]-phibar[j,-7]) # cond. mean phi8i|rest
      v<-(sigmaphi[7,7]-sigmaphi[7,-7]%*%solve(sigmaphi[-7,-7])%*%sigmaphi[-7,7])/mi[j]    # cond. var
      pnew<-dnorm(phi7new[j],m, sqrt(v), log=T)
      pold<-dnorm(phi7[j],m, sqrt(v), log=T)
      ratio<-(lnew[j]+pnew)-(lold[j]+pold)
      if (log(runif(1))<ratio) {
        phi7[j]<-phi7new[j]
        Accu3[j]<-Accu3[j]+1  # keep track of acceptance of phi for each county
      }
    }
    
    # Center phi7 
    phi7<-phi7-mean(phi7)
    Phi7<-rep(phi7,nis)
    phimat[,7]<-phi7
    
    #---------#
    # BB phi8 #
    #---------#
    phi8new<-phi8+s8*rt(n, 3)                   # Proposal distribution
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+rep(phi8new,nis)*int
    mu<-1/(1+exp(-eta))
    lnew<-rep(0,n)                              # account for empty counties - their likelihood contribution is 0
    lnew[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    eta<-X%*%beta+Phi5+Phi6*t+Phi7*nhb+Phi8*int
    mu<-1/(1+exp(-eta))
    lold<-rep(0,n)                              # account for empty counties - their likelihood contribution is 0
    lold[nis1>0]<-tapply(dbetabinom(y[y1==1],ntrial,mu[y1==1],rho,log=T),id[y1==1],sum)
    
    for (j in 1:n){
      phibar[j,8]<-mean(phi8[which(A[j,]==1)])  # update phi8bar[j] before updating phi8[j]
      m<-phibar[j,8]+sigmaphi[8,-8]%*%solve(sigmaphi[-8,-8])%*%(phimat[j,-8]-phibar[j,-8])  # cond. mean phi8i|rest
      v<-(sigmaphi[8,8]-sigmaphi[8,-8]%*%solve(sigmaphi[-8,-8])%*%sigmaphi[-8,8])/mi[j]     # cond. var
      pnew<-dnorm(phi8new[j],m, sqrt(v), log=T)
      pold<-dnorm(phi8[j],m, sqrt(v), log=T)
      ratio<-(lnew[j]+pnew)-(lold[j]+pold)
      if (log(runif(1))<ratio) {
        phi8[j]<-phi8new[j]
        Accu4[j]<-Accu4[j]+1                   # keep track of acceptance of phi for each county
      }
    }
    
    # Center phi8
    phi8<-phi8-mean(phi8)
    Phi8<-rep(phi8,nis)
    phimat[,8]<-phi8
    
    #-----------------#
    # Update sigmaphi #
    #-----------------#
    phimat<-cbind(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8)
    sigmaphi<-riwish(c0+n-1,d0+t(phimat)%*%Q%*%phimat)
    
    #---------------#
    # Store Results #
    #---------------#
    if (i> burn & i%%thin==0) {
      j<-(i-burn)/thin
      Alpha[j,]<-alpha
      Beta[j,]<-beta
      Rho[j]<-rho
      S2Phi[j,]<-c(sigmaphi)
      Phi1s[j,]<-phi1
      Phi2s[j,]<-phi2
      Phi3s[j,]<-phi3
      Phi4s[j,]<-phi4
      Phi5s[j,]<-phi5
      Phi6s[j,]<-phi6
      Phi7s[j,]<-phi7
      Phi8s[j,]<-phi8
    }
    if (i%%10==0)  {
      print(i)
      print(alpha)
      print(beta)
      print(rho)
      print(round(sigmaphi,3))
    }
  }
  (time.tol<-proc.time()-time.start)
  
  #---------#
  # Counter #
  #---------#
  A1/nsim
  A2/nsim
  Accu1/nsim
  Accu2/nsim
  Accu3/nsim
  Accu4/nsim
  
  #---------#
  # Results #
  #---------#
  malpha<-colMeans(Alpha)
  qalpha<-apply(Alpha,2,quantile,c(0.025,0.975))
  
  mbeta<-colMeans(Beta)
  qbeta<-apply(Beta,2,quantile,c(0.025,0.975))
  
  mrho<-mean(Rho)
  qrho<-quantile(Rho,c(0.025,0.975))
  
  ms2phi<-colMeans(S2Phi)
  qs2phi<-apply(S2Phi,2,quantile,c(0.025,0.975))
  
  true.alpha
  malpha
  qalpha
  true.beta
  mbeta
  qbeta
  mrho
  qrho
  cov
  ms2phi
  qs2phi
  
  
  #------#
  # WAIC #
  #------#
  L<-matrix(0,lastit,N)  # Select the last 5000 iters
  for (j in 1:lastit){
    # likelihood function for ZIBB
    eta1<-X%*%Alpha[j,]+rep(Phi1s[j,],nis)+rep(Phi2s[j,],nis)*t+rep(Phi3s[j,],nis)*nhb+rep(Phi4s[j,],nis)*int
    mu1<-1/(1+exp(-eta1))
    
    eta2<-X%*%Beta[j,]+rep(Phi5s[j,],nis)+rep(Phi6s[j,],nis)*t+rep(Phi7s[j,],nis)*nhb+rep(Phi8s[j,],nis)*int
    mu2<-1/(1+exp(-eta2))
    
    L[j,]<-dZIBB(y,mu=mu2,sigma=Rho[j]/(1-Rho[j]),nu=1-mu1,bd=ntrial)  # dZIBB() 
    
    if (j%%100==0) print(j)
  }
  lhat<-colMeans(L)
  lppd<-sum(log(lhat))
  pwaic<-sum(apply(log(L),2,var))
  waic<- -2*(lppd-pwaic)
  cat("WAIC:",waic)  # WAIC: 250584.7
  
  #-------------#
  # Trace Plots #
  #-------------#
  plot(1:lastit,Alpha[,1],type="l",col="darkgreen",xlab="Iteration",ylab=expression(alpha[1]))
  abline(h=malpha[1],col="blue4")
  plot(1:lastit,Alpha[,2],type="l",col="darkgreen",xlab="Iteration",ylab=expression(alpha[2]))
  abline(h=malpha[2],col="blue4")
  plot(1:lastit,Alpha[,3],type="l",col="darkgreen",xlab="Iteration",ylab=expression(alpha[3]))
  abline(h=malpha[3],col="blue4")
  plot(1:lastit,Alpha[,4],type="l",col="darkgreen",xlab="Iteration",ylab=expression(alpha[4]))
  abline(h=malpha[4],col="blue4")
  
  
  plot(1:lastit,Beta[,1],type="l",col="darkgreen",xlab="Iteration",ylab=expression(beta[1]))
  abline(h=mbeta[1],col="blue4")
  plot(1:lastit,Beta[,2],type="l",col="darkgreen",xlab="Iteration",ylab=expression(beta[2]))
  abline(h=mbeta[2],col="blue4")
  plot(1:lastit,Beta[,3],type="l",col="darkgreen",xlab="Iteration",ylab=expression(beta[3]))
  abline(h=mbeta[3],col="blue4")
  plot(1:lastit,Beta[,4],type="l",col="darkgreen",xlab="Iteration",ylab=expression(beta[4]))
  abline(h=mbeta[4],col="blue4")
  
  plot(1:lastit,Rho,type="l",col="darkgreen",xlab="Iteration",ylab=expression(rho))
  abline(h=mrho,col="blue4")
  
  plot(1:lastit,S2Phi[,1],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi1]))
  abline(h=ms2phi[1],col="blue4")
  plot(1:lastit,S2Phi[,10],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi2]))
  abline(h=ms2phi[10],col="blue4")
  plot(1:lastit,S2Phi[,19],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi3]))
  abline(h=ms2phi[19],col="blue4")
  plot(1:lastit,S2Phi[,28],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi4]))
  abline(h=ms2phi[28],col="blue4")
  plot(1:lastit,S2Phi[,37],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi5]))
  abline(h=ms2phi[37],col="blue4")
  plot(1:lastit,S2Phi[,46],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi6]))
  abline(h=ms2phi[46],col="blue4")
  plot(1:lastit,S2Phi[,55],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi7]))
  abline(h=ms2phi[55],col="blue4")
  plot(1:lastit,S2Phi[,64],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi8]))
  abline(h=ms2phi[64],col="blue4")
  
  plot(1:lastit,S2Phi[,5],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi12]))
  abline(h=ms2phi[5],col="blue4")
  plot(1:lastit,S2Phi[,9],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi13]))
  abline(h=ms2phi[9],col="blue4")
  plot(1:lastit,S2Phi[,10],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi23]))
  abline(h=ms2phi[10],col="blue4")
  plot(1:lastit,S2Phi[,13],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi14]))
  abline(h=ms2phi[13],col="blue4")
  plot(1:lastit,S2Phi[,14],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi24]))
  abline(h=ms2phi[14],col="blue4")
  plot(1:lastit,S2Phi[,15],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi34]))
  abline(h=ms2phi[15],col="blue4")
  
  
  # Radom selected county-specific effects
  plot(1:lastit,Phi3s[,20],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi3s]))
  plot(1:lastit,Phi4s[,15],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi4s]))
  
  plot(1:lastit,Phi1s[,5],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi3s]))
  plot(1:lastit,Phi2s[,20],type="l",col="darkgreen",xlab="Iteration",ylab=expression(Sigma[phi3s]))
  
  #------------------------------------------------------------------------------#
  
  #-------------------#
  # Save MCMC Samples #
  #-------------------#
  
  samples<-list(Alpha=Alpha,
                Beta=Beta,
                Rho=Rho,
                S2Phi=S2Phi,
                Phi1s=Phi1s,
                Phi2s=Phi2s,
                Phi3s=Phi3s,
                Phi4s=Phi4s,
                Phi5s=Phi5s,
                Phi6s=Phi6s,
                Phi7s=Phi7s,
                Phi8s=Phi8s)
  
  # dir.sav<-"..\\"
  # save(samples,file=paste(dir.sav,"Simulation-ZIBBSVC.Rda",sep=""))  # Save MCMC samples file
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
  
  # Binary Random Intercept (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi1"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[1*i*1], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Intercept (Estiamted)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi1p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[1*i*1], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Slope (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi2"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[1*i*2], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Slope (ESTIAMTED)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi2p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[1*i*2], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Race (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi3"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[1*i*3], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Race (ESTIAMTED) 
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi3p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[1*i*3], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Interaction (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi4"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[1*i*4], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # Binary Random Interaction (ESTIMATED)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi4p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[1*i*4], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Intercept (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi5"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[2*i*1], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Intercept (ESTIMATED)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi5p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[2*i*1], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Slope (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi6"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[2*i*2], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Slope (ESTIMATED)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi6p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[2*i*2], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Race (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi7"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[2*i*3], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Race (ESTIAMTED)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi7p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[2*i*3], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Interaction (TRUE)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi8"), midpoint = c(NA), title =  expression(paste("Simulated ", phi[2*i*4], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  # BB Random Interaction (ESTIAMTED)
  tm_shape(sc_mcmf_data)+
    tm_fill(c("phi8p"), midpoint = c(NA), title =  expression(paste("Estimated ", phi[2*i*4], sep = " ")), palette = pal, style = "quantile")+
    tm_layout(title = "",
              title.size = 4.0,
              title.position = c("right", "bottom"),
              title.fontface = "bold",
              legend.outside.size=0.25,
              legend.text.size = 1.5,
              legend.title.fontface = "bold",
              main.title.fontface = "bold",
              main.title.position = "center",
              legend.title.size = 2,
              panel.label.size = 2.5,
              legend.format = c(digits=2))+
    tm_borders(alpha = 0.3, lwd = 1)
  
  
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
  
  ggplot(dcounty,aes(x=t,y=mu,group=gp,col=gp))+        # all 46 county-specific trends
    geom_line()+
    geom_ribbon(aes(ymin = lb, ymax = ub,col=gp,fill=gp,group=gp),linetype="twodash",alpha=0.4,show.legend = F)+
    facet_wrap(~countyid)    
  
    
  
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
  
