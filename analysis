##=============================================================
# Function: analysis
#
# Description:
#   Estimates the model parameters for clustered interval-censored
#     competing-risks data in the presence of informative cluster size (ICS)
#     using B-spline-based sieve estimation.
#
# Model:
#   Generalized odds-rate transformation model
#
# Input:
#   simulation --- a data frame containing:
#     identi    -  the cluster identifier for the subject  
#     id        -  an arbitrary subject index within the cluster; its ordering is not important
#     sizenum   -  the total number of subjects in the corresponding cluster
#     v         -  the left endpoint of the observation interval
#     u         -  the right endpoint of the observation interval
#                   (Inf is allowed for right-censored observations)
#     J         -  the event indicator
#     X1,...,Xp -  covariates included in the model
#
#   alpha1      -  the transformation parameter for cause 1
#   alpha2      -  the transformation parameter for cause 2
#
# Covariates:
#   The function allows an arbitrary number of covariates.
#   However, we recommend using at most five covariates to avoid
#     exponential growth in the number of boundary combinations.
#
# Output:
#   The convergence status of the optimization is first printed to the console:
#     "Convergence successful" or "Convergence failed".
#
#   The function then returns a numeric vector containing, in order:
#    (1) regression coefficient estimates for cause 1
#    (2) regression coefficient estimates for cause 2
#    (3) estimated baseline function values for cause 1
#    (4) estimated baseline function values for cause 2
#
#
# Required package:
#   alabama
#   splines
#
#==============================================================


analysis<-function(simulation,alpha1,alpha2){
  
  data<-simulation
  n<-length(data$identi)
 
  variable<-grep("^X[0-9]+$",names(data),value=TRUE)
  
  order<-order(as.integer(sub("^X","",variable)))
  variable<-variable[order]
 
  X<-as.matrix(data[,variable,drop=FALSE])
  qc<-ncol(X)
  
  weight<-1/(data$sizenum)
  
  tao_1<-min(data$v)
  datau<-data$u[is.finite(data$u)]
  tao_2<-max(datau)
  t<-c(data$v,datau)
  splinenum<-floor(length(unique(t))^(1/3))
#  splinenum<-floor(length(unique(t))^(1/5))
  tao<-c(tao_1,tao_2)
   
  knots<-unique(quantile(t[t<max(t)&t>min(t)],seq(0,1,by=1/(splinenum+1)))[2:(splinenum+1)])
  baseBfunction<-splines::bs(t,knots=knots,degree=3,intercept=TRUE,Boundary.knots=c(min(t),max(t)))     
 
  BfuncV<-predict(baseBfunction,data$v)
   
  delta<-data$J
  TU<-data$u
  TU[delta==0]<-max(t)
  BfuncU<-predict(baseBfunction,TU)
  
  ns <-dim(BfuncU)[2]
  
  wtheta1<-rep(NA,qc)
  wtheta2<-rep(NA,qc)
  wbeta1<-rep(NA,ns)
  wbeta2<-rep(NA,ns)
  wparameter<-c(wtheta1,wtheta2,wbeta1,wbeta2)
  wparameter<-as.numeric(wparameter)
  
  Tv<-data$v
  
  d1<-(delta == 1 & Tv > 0)
  d2<-(delta == 2 & Tv > 0)
  d1_1<-(delta == 1 & Tv == 0)
  d2_1<-(delta == 2 & Tv == 0)
  d<-(d1 + d1_1 + d2 + d2_1)
  
  ## Generate original beta
  origin<-seq(from = 0.0001, to = 0.875, by = ((0.875 - 0.0001) / (dim(baseBfunction)[2] - 1)))
  origin<-log(origin ^ 3)
  theta<-c(origin,origin)
  originpara<-c(rep(0,2*qc),theta)
  
  loglikelihood<-function(para){
    #辅助变量x<-originpara
    x<-para
    ######attention!########
    wtheta1<-x[1:qc]
    wtheta2<-x[(qc+1):(2*qc)]
    wbeta1<-x[(2*qc+1):(2*qc+ns)]
    wbeta2<-x[(2*qc+ns+1):(2*qc+2*ns)]
    
    BS1u <- BfuncU %*% wbeta1
    BS1v <- BfuncV %*% wbeta1
    BS2u <- BfuncU %*% wbeta2
    BS2v <- BfuncV %*% wbeta2
    bz_1 <- X %*% wtheta1
    bz_2 <- X %*% wtheta2
    
    if(alpha1>0){
      ci1v <- 1 - (1 +alpha1*exp(BS1v + bz_1))^(-1/alpha1)
      ci1u <- 1 - (1 +alpha1*exp(BS1u + bz_1))^(-1/alpha1)  
    }else if(alpha1==0){
      ci1v <- 1-exp(-exp(BS1v + bz_1))
      ci1u <- 1-exp(-exp(BS1u + bz_1))
    }
    
    if(alpha2>0){
      ci2v <- 1 - (1 +alpha2*exp(BS2v + bz_2))^(-1/alpha2)
      ci2u <- 1 - (1 +alpha2*exp(BS2u + bz_2))^(-1/alpha2)  
    }else if(alpha2==0){
      ci2v <- 1-exp(-exp(BS2v + bz_2))
      ci2u <- 1-exp(-exp(BS2u + bz_2))
    }
    
    ci1u[ci1u==ci1v&delta==0]<-ci1u[ci1u==ci1v&delta==0]+0.001
    ci2u[ci2u==ci2v&delta==0]<-ci2u[ci2u==ci2v&delta==0]+0.001
    
    ill <- (d1_1 * log(ci1u) + d2_1 * log(ci2u) +
              d1 * log(ci1u - ci1v) + d2 * log(ci2u - ci2v) +
              (1 - d) * log(1 - (ci1v + ci2v)))*weight
    
    nll <- -sum(ill)
    nll
  }
  
  Gradfunction<-function(para){
    x<-para
    ######attention!########
    wtheta1<-x[1:qc]
    wtheta2<-x[(qc+1):(2*qc)]
    wbeta1<-x[(2*qc+1):(2*qc+ns)]
    wbeta2<-x[(2*qc+ns+1):(2*qc+2*ns)]
    
    BS1u <- as.vector(BfuncU %*% wbeta1)
    BS1v <- as.vector(BfuncV %*% wbeta1)
    BS2u <- as.vector(BfuncU %*% wbeta2)
    BS2v <- as.vector(BfuncV %*% wbeta2)
    bz_1 <- as.vector(X %*% wtheta1)
    bz_2 <- as.vector(X %*% wtheta2)
    
    if(alpha1>0){
      ci1v <- 1-(1+alpha1*exp(BS1v + bz_1))^(-1/alpha1)
      ci1u <- 1-(1+alpha1*exp(BS1u + bz_1))^(-1/alpha1)  
    }else if(alpha1==0){
      ci1v <- 1-exp(-exp(BS1v+bz_1))
      ci1u <- 1-exp(-exp(BS1u+bz_1))
    }
    
    if(alpha2>0){
      ci2v <- 1-(1+alpha2*exp(BS2v + bz_2))^(-1/alpha2)
      ci2u <- 1-(1+alpha2*exp(BS2u + bz_2))^(-1/alpha2)  
    }else if(alpha2==0){
      ci2v <- 1-exp(-exp(BS2v+bz_2))
      ci2u <- 1-exp(-exp(BS2u+bz_2))
    } 
    ci1u[ci1u==ci1v&delta==0]<-ci1u[ci1u==ci1v&delta==0]+0.001
    ci2u[ci2u==ci2v&delta==0]<-ci2u[ci2u==ci2v&delta==0]+0.001
    
    zero <- matrix(0,nrow=n,ncol=qc)
    dB1u <- BfuncU
    dB1u <- cbind(dB1u,matrix(0,nrow=n,ncol=ns))
    dB1u <- cbind(X,zero,dB1u)
    dB2u <- BfuncU
    dB2u <- cbind(matrix(0,nrow=n,ncol=ns),dB2u)
    dB2u <- cbind(zero,X,dB2u)
    
    dB1v <- BfuncV
    dB1v <- cbind(dB1v,matrix(0,nrow=n,ncol=ns))
    dB1v <- cbind(X,zero,dB1v)
    dB2v <- BfuncV
    dB2v <- cbind(matrix(0,nrow=n,ncol=ns),dB2v)
    dB2v <- cbind(zero,X,dB2v)
    
    if(alpha1>0){
      dci1v <- (1+alpha1*exp(BS1v+bz_1))^(-(1/alpha1)-1)*exp(BS1v+bz_1)
      dci1u <- (1+alpha1*exp(BS1u+bz_1))^(-(1/alpha1)-1)*exp(BS1u+bz_1)  
    }else if(alpha1==0){
      dci1v <- exp(-exp(BS1v + bz_1)) * exp(BS1v + bz_1)
      dci1u <- exp(-exp(BS1u + bz_1)) * exp(BS1u + bz_1)
      
    }
    
    if(alpha2>0){
      dci2v <- (1+alpha2*exp(BS2v+bz_2))^(-(1/alpha2)-1)*exp(BS2v+bz_2)
      dci2u <- (1+alpha2*exp(BS2u+bz_2))^(-(1/alpha2)-1)*exp(BS2u+bz_2)
      
    }else if(alpha2==0){
      dci2v <- exp(-exp(BS2v + bz_2)) * exp(BS2v + bz_2)
      dci2u <- exp(-exp(BS2u + bz_2)) * exp(BS2u + bz_2)
    }
    
    gradlikeli <- ((d1_1 / ci1u) * (dci1u*dB1u) +
                     (d2_1 / ci2u) * (dci2u * dB2u) +
                     (d1 / (ci1u - ci1v)) * (dci1u* dB1u - dci1v* dB1v) +
                     (d2 / (ci2u - ci2v)) * (dci2u* dB2u - dci2v* dB2v) +
                     ((1 - d) / (1 - ci1v - ci2v)) * (-dci1v* dB1v - dci2v * dB2v))*weight
   
    sumgradlikeli<- -colSums(gradlikeli) 
    sumgradlikeli
  }
  
  if(qc==1){
    comb<-matrix(c(min(X), max(X)), ncol = 1)
  }else{
    mM<-rbind(apply(X, 2, min), apply(X, 2, max))
    comp <- function(x){
      mM[,x]
    }
    comb <- expand.grid(lapply(1:ncol(X), comp))
  }
  
  comb<-as.matrix(comb)
  
  hin<-function(para){
    #para<-originpara
    x<-para
    
    wtheta1<-x[1:qc]
    wtheta2<-x[(qc+1):(2*qc)]
    wbeta1<-x[(2*qc+1):(2*qc+ns)]
    wbeta2<-x[(2*qc+ns+1):(2*qc+2*ns)]
   
    wminus1<-diff(wbeta1)
    wminus2<-diff(wbeta2)
   
    #CIF
    
    cif1 <- function(i,wtheta1,wbeta1){
      if(alpha1 > 0){
        (1 + alpha1* exp(sum(comb[i,]*wtheta1)+wbeta1[length(wbeta1)]))^(-1 /alpha1)
      } else if(alpha1== 0){
        exp(-exp(sum(comb[i,]*wtheta1)+wbeta1[length(wbeta1)]))
      }
    }
    cif2 <- function(i,wtheta2,wbeta2){
      if(alpha2 > 0){
        (1 + alpha2 * exp(sum(comb[i,]*wtheta2)+wbeta2[length(wbeta2)]))^(-1 /alpha2 )
      } else if(alpha2 == 0){
        exp(-exp(sum(comb[i,]*wtheta2)+wbeta2[length(wbeta2)]))
      }
    }
    addcons<-rep(NA,length(comb[,1]))
    for(j in 1:length(comb[,1])){
      addcons[j]<-cif1(j,wtheta1,wbeta1)+cif2(j,wtheta2,wbeta2)-1
    }
    ui<-c(wminus1,wminus2,addcons)
    ui-0.00000001
  }
  
  
  hin_grad<-function(para){
    #para<-originpara
    x<-para
    
    wtheta1<-x[1:qc]
    wtheta2<-x[(qc+1):(2*qc)]
    wbeta1<-x[(2*qc+1):(2*qc+ns)]
    wbeta2<-x[(2*qc+ns+1):(2*qc+2*ns)]
    wminus1<-diff(wbeta1)
    wminus2<-diff(wbeta2)
    
    zeros1<-matrix(0,nrow=length(wminus1)+length(wminus2),ncol=length(wtheta1)+length(wtheta2))
    zeros2<-matrix(0,nrow=length(wminus1),ncol=length(wbeta1))
    
    
    s<-length(wbeta1)
     
    u<-matrix(0, nrow=s-1, ncol=s)
     
    for (i in 1:(s-1)){
      u[i,i]<- -1
      u[i,i+1]<- 1
    }
    transi1<-cbind(u,zeros2)
    transi2<-cbind(zeros2,u)
    transi3<-rbind(transi1,transi2)
    ui<-cbind(zeros1,transi3)
    
    
    dcif1 <- function(i,wtheta1,wbeta1){
      if(alpha1 > 0){
        -(1 + alpha1* exp(sum(comb[i,]*wtheta1)+wbeta1[length(wbeta1)]))^(-(1 /alpha1) - 1) * exp(sum(comb[i,]*wtheta1)+wbeta1[length(wbeta1)])
      } else if(alpha1 == 0){
        -exp(-exp(sum(comb[i,]*wtheta1)+wbeta1[length(wbeta1)]))*exp(sum(comb[i,]*wtheta1)+wbeta1[length(wbeta1)])
      }
    }
    
    dcif2 <- function(i,wtheta2,wbeta2){
      if(alpha2 > 0){
        -(1 + alpha2* exp(sum(comb[i,]*wtheta2)+wbeta2[length(wbeta2)]))^(-(1 / alpha2) - 1) * exp(sum(comb[i,]*wtheta2)+wbeta2[length(wbeta2)])
      } else if(alpha2== 0){
        -exp(-exp(sum(comb[i,]*wtheta2)+wbeta2[length(wbeta2)]))*exp(sum(comb[i,]*wtheta2)+wbeta2[length(wbeta2)])
      }
    }
    
    addcons<-rep(0,length(wparameter))
    
    for(k in 1:length(comb[,1])){
      
      for(j in 1:qc){
        addcons[j]<-dcif1(k,wtheta1,wbeta1)*comb[k,j]
        addcons[qc+j]<-dcif2(k,wtheta2,wbeta2)*comb[k,j]
      }
      addcons[2*qc+ns]<-dcif1(k,wtheta1,wbeta1)
      addcons[2*qc+2*ns]<-dcif2(k,wtheta2,wbeta2)
      ui<-rbind(ui,addcons)
    
    }
    unname(ui)
  }
  
  
  heq<-function(para){
   
    x<-para
    
    wtheta1<-x[1:qc]
    wtheta2<-x[(qc+1):(2*qc)]
    wbeta1<-x[(2*qc+1):(2*qc+ns)]
    wbeta2<-x[(2*qc+ns+1):(2*qc+2*ns)]
    wminus1<-diff(wbeta1)
    wminus2<-diff(wbeta2)
   
    wminus1<-rep(0,ns-1)
    wminus2<-rep(0,ns-1)
     
    
    cif1 <- function(i,wtheta1,wbeta1){
      if(alpha1 > 0){
        (1 + alpha1 * exp(sum(comb[i,]*wtheta1)+wbeta1[1]))^(-1 / alpha1 )
      } else if(alpha1== 0){
        exp(-exp(sum(comb[i,]*wtheta1)+wbeta1[1]))
      }
    }
    cif2 <- function(i,wtheta2,wbeta2){
      if(alpha2 > 0){
        (1 + alpha2 * exp(sum(comb[i,]*wtheta2)+wbeta2[1]))^(-1 / alpha2)
      } else if(alpha2 == 0){
        exp(-exp(sum(comb[i,]*wtheta2)+wbeta2[1]))
        
      }
    }
    addcons<-rep(NA,length(comb[,1]))
    for(j in 1:length(comb[,1])){
      addcons[j]<-cif1(j,wtheta1,wbeta1)+cif2(j,wtheta2,wbeta2)-2
    }
    ui<-c(wminus1,wminus2,addcons)
    ui    
    
  }
  
  
  heq_grad<-function(para){
    x<-para
    
    wtheta1<-x[1:qc]
    wtheta2<-x[(qc+1):(2*qc)]
    wbeta1<-x[(2*qc+1):(2*qc+ns)]
    wbeta2<-x[(2*qc+ns+1):(2*qc+2*ns)]
    
    ui<-matrix(0,nrow=2*(ns-1),ncol=length(x))
    dcif1 <- function(i,wtheta1,wbeta1){
      if(alpha1 > 0){
        (1 +alpha1* exp(sum(comb[i,]*wtheta1)+wbeta1[1]))^(-(1 /alpha1) - 1) * exp(sum(comb[i,]*wtheta1)+wbeta1[1])
      } else if(alpha1== 0){
        exp(-exp(sum(comb[i,]*wtheta1)+wbeta1[1])) * exp(sum(comb[i,]*wtheta1)+wbeta1[1])
        
      }
    }
    dcif2 <- function(i,wtheta2,wbeta2){
      if(alpha2 > 0){
        (1 + alpha2* exp(sum(comb[i,]*wtheta2)+wbeta2[1]))^(-(1 /alpha2) - 1) * exp(sum(comb[i,]*wtheta2)+wbeta2[1])
      } else if(alpha2 == 0){
        exp(-exp(sum(comb[i,]*wtheta2)+wbeta2[1])) * exp(sum(comb[i,]*wtheta2)+wbeta2[1])
      }
    }
    
    addcons<-rep(0,length(wparameter))
    for(k in 1:length(comb[,1])){
      for(j in 1:qc){
        addcons[j]<- dcif1(k,wtheta1,wbeta1)*comb[k,j]
        addcons[qc+j]<- dcif2(k,wtheta2,wbeta2)*comb[k,j]
      }
      addcons[2*qc+1]<- dcif1(k,wtheta1,wbeta1)
      addcons[2*qc+ns+1]<- dcif2(k,wtheta2,wbeta2)
     
      ui<-rbind(ui,addcons)
      
    }
    unname(ui)
  }
  
  
  others <- try(alabama::constrOptim.nl(par = originpara,
                                        fn = loglikelihood,
                                        gr = Gradfunction,
                                        hin = hin,
                                        hin.jac = hin_grad,
                                        heq = heq,
                                        heq.jac = heq_grad,
                                        control.optim = list(maxit = 2000),
                                        control.outer = list(trace = FALSE)), silent = TRUE)
  
  
  if (others$convergence == 0) {
    print("Convergence successful")
  } else {
    print("Convergence failed")
  }
  
  coef <- others$par[1:(2*qc)]
  Bfuncpre <- predict(baseBfunction, 1:250/50)
  estivalue1 <- Bfuncpre %*% others$par[(2*qc+1):(2*qc+ns)]
  estivalue2 <- Bfuncpre %*% others$par[(2*qc+ns+1):(2*qc+2*ns)]
  every <- c(coef,exp(estivalue1), exp(estivalue2))
  
  return(every)
  
}


