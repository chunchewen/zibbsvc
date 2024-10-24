  # Zero-Inflated Beta-Binomial Spatially Varying Coefficients (ZIBB-SVC) Model
  # Spatially varying space x time in both component
  # Simulation Analysis: Generate Simulated Data and Run MCMC 
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

  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  
