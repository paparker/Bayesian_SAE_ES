library(Matrix)

mod1Fit <- function(Y, X, Var, sig2B=100, iter=100, burn=50){
  n <- length(Y)
  p <- ncol(X)
  
  sig2RE <- 1
  eta <- rep(0,n)
  beta <- rep(1,p)
  etaOut <- matrix(NA, nrow=iter, ncol=n)
  betaOut <- matrix(NA, nrow=iter, ncol=p)
  
  preds <- matrix(NA, nrow=n, ncol=iter)
  
  Binv <- Diagonal(p, 1/sig2B)
  precBeta <- t(X)%*%Diagonal(x=1/Var)%*%X+Binv
  Ub <- chol(forceSymmetric(precBeta))
  
  pb=txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## Fixed effects
    meanBeta <- t(X)%*%Diagonal(x=1/Var)%*%(Y-eta)
    b <- rnorm(nrow(precBeta))
    beta <- betaOut[i,] <- backsolve(Ub, backsolve(Ub, meanBeta, transpose=T) + b)
    
    ## Random effects
    precEta <- Diagonal(x=1/Var+(1/sig2RE))
    U <- chol(forceSymmetric(precEta))
    meanEta <- Diagonal(x=1/Var)%*%(Y-X%*%beta)
    b <- rnorm(nrow(precEta))
    eta <- etaOut[i,] <- backsolve(U, backsolve(U, meanEta, transpose=T) + b)
    
    ## RE variance
    sig2RE <- 1/rgamma(1, shape=(n-1)/2, (0.5 + 0.5*sum(eta^2)))
    
    ## Predictions
    preds[,i] <- X%*%beta + eta
    
    setTxtProgressBar(pb, i)
  }
  return(list(Preds=preds[,-c(1:burn)], Beta=betaOut[-c(1:burn),], Eta=etaOut[-c(1:burn),]))
}


mod2Fit <- function(Y, X, Psi1, Psi2, Var, sig2B=100, iter=100, burn=50){
  n <- length(Y)
  p <- ncol(X)
  r1 <- ncol(Psi1)
  r2 <- ncol(Psi2)
  C <- cbind(X,Psi1,Psi2)
  
  
  sig2RE1 <- 1
  sig2RE2 <- 1
  eta1 <- rep(0,r1)
  eta2 <- rep(0,r2)
  beta <- rep(1,p)
  eta1Out <- matrix(NA, nrow=iter, ncol=r1)
  eta2Out <- matrix(NA, nrow=iter, ncol=r2)
  betaOut <- matrix(NA, nrow=iter, ncol=p)
  
  preds <- matrix(NA, nrow=n, ncol=iter)
  
  Binv <- Diagonal(p, 1/sig2B)
  precBeta <- t(X)%*%Diagonal(x=1/Var)%*%X+Binv
  Ub <- chol(forceSymmetric(precBeta))
  
  pb=txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    Einv1 <- Diagonal(r1, 1/sig2RE1)
    Einv2 <- Diagonal(r2, 1/sig2RE2)
    precAll <- t(C)%*%Diagonal(x=1/Var)%*%C+bdiag(Binv,Einv1,Einv2) 
    U <- chol(forceSymmetric(precAll))
    meanAll <- t(C)%*%Diagonal(x=1/Var)%*%Y
    b <- rnorm(nrow(precAll))
    tmp <- backsolve(U, backsolve(U, meanAll, transpose=T) + b)
    beta <- betaOut[i,] <- tmp[1:p]
    eta1 <- eta1Out[i,] <- tmp[(p+1):(p+r1)]
    eta2 <- eta2Out[i,] <- tmp[-c(1:(p+r1))]
      
    
    ## RE variance
    sig2RE1 <- 1/rgamma(1, shape=(r1-1)/2, (0.5 + 0.5*sum(eta1^2)))
    sig2RE2 <- 1/rgamma(1, shape=(r2-1)/2, (0.5 + 0.5*sum(eta2^2)))
    
    ## Predictions
    preds[,i] <- X%*%beta + Psi1%*%eta1 + Psi2%*%eta2
    
    setTxtProgressBar(pb, i)
  }
  return(list(Preds=preds[,-c(1:burn)], 
              Beta=betaOut[-c(1:burn),], 
              Eta1=eta1Out[-c(1:burn),], 
              Eta2=eta2Out[-c(1:burn),]))
}

mod3Fit <- function(Y, X, Psi1, Psi2, wgt, predX, predPsi1, predPsi2, predGrp, sig2B=100, iter=100, burn=50){
  n <- length(Y)
  nP <- nrow(predX)
  p <- ncol(X)
  r1 <- ncol(Psi1)
  r2 <- ncol(Psi2)
  C <- cbind(X,Psi1,Psi2)
  
  sig2y <- 1
  sig2RE1 <- 1
  sig2RE2 <- 1
  eta1 <- rep(0,r1)
  eta2 <- rep(0,r2)
  beta <- rep(1,p)
  eta1Out <- matrix(NA, nrow=iter, ncol=r1)
  eta2Out <- matrix(NA, nrow=iter, ncol=r2)
  betaOut <- matrix(NA, nrow=iter, ncol=p)
  sig2yOut <- rep(NA,iter)
  
  preds <- matrix(NA, nrow=length(predGrp), ncol=iter-burn)
  
  Binv <- Diagonal(p, 1/sig2B)
  
  pb=txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    Einv1 <- Diagonal(r1, 1/sig2RE1)
    Einv2 <- Diagonal(r2, 1/sig2RE2)
    precAll <- t(C)%*%Diagonal(x=wgt/sig2y)%*%C+bdiag(Binv,Einv1,Einv2) 
    U <- chol(forceSymmetric(precAll))
    meanAll <- t(C)%*%Diagonal(x=wgt/sig2y)%*%Y
    b <- rnorm(nrow(precAll))
    tmp <- backsolve(U, backsolve(U, meanAll, transpose=T) + b)
    beta <- betaOut[i,] <- tmp[1:p]
    eta1 <- eta1Out[i,] <- tmp[(p+1):(p+r1)]
    eta2 <- eta2Out[i,] <- tmp[-c(1:(p+r1))]
    
    
    ## RE variance
    sig2y <- sig2yOut[i] <- 1/rgamma(1, shape=(n+1)/2, 0.5 + 0.5*sum((wgt*(Y-X%*%beta-Psi1%*%eta1-Psi2%*%eta2))^2))
    sig2RE1 <- 1/rgamma(1, shape=(r1+1)/2, (0.5 + 0.5*sum(eta1^2)))
    sig2RE2 <- 1/rgamma(1, shape=(r2+1)/2, (0.5 + 0.5*sum(eta2^2)))
    
    ## Predictions
    if(i>burn){
      preds[,i-burn] <- rnorm(nP, mean=as.numeric(predX%*%beta + predPsi1%*%eta1 + predPsi2%*%eta2), sd=sqrt(sig2y))
    }
    #tmp <- exp(rnorm(nP, mean=as.numeric(predX%*%beta + predPsi1%*%eta1 + predPsi2%*%eta2), sd=sqrt(sig2y)))
    #tmpDF <- data.frame(Grp=predGrp, tmp) %>% group_by(Grp) %>% summarize(Est=mean(tmp), .groups='drop') 
    #preds[,i] <- tmpDF$Est
    setTxtProgressBar(pb, i)
  }
  predDF <- data.frame(Grp=predGrp, preds) %>% group_by(Grp) %>% summarize_all(mean)
  return(list(Preds=predDF, 
              Beta=betaOut[-c(1:burn),], 
              Eta1=eta1Out[-c(1:burn),], 
              Eta2=eta2Out[-c(1:burn),],
              sig2y=sig2yOut))
}
