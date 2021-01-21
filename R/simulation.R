library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mase)
library(ggthemes)
library(sampling)
library(sae)
source('R/saeFuncs.R')

## Simulation settings
Nsim <- 100
ss <- 20000

## Ready survey data
pums <- read_csv('Data/pums.csv')
pums <- pums %>% dplyr::select(SECTOR, FIPST, TABWGT, EMP=EMPLOYMENT_NOISY,
                               PAYROLL=PAYROLL_NOISY, RECEIPTS=RECEIPTS_NOISY,
                               ECOM=WEBSITE, HOMEBASED)
pums <- pums %>% mutate(ECOM=as.numeric(ECOM==1))
pums$ECOM[is.na(pums$ECOM)] <- 0

truth <- pums %>% group_by(SECTOR, FIPST) %>% summarize(Mean=mean(RECEIPTS), Size=n())
truth$Grp <- 1:nrow(truth)
pums <- pums %>% left_join(truth%>%dplyr::select(SECTOR,FIPST,Grp), by=c("SECTOR","FIPST"))
predXT <- pums %>% group_by(SECTOR, FIPST) %>% summarize(EMP=mean(EMP), ECOM=mean(ECOM))

HTests <- UWests <- mod1ests <- mod2ests <- mod3ests <- matrix(NA, nrow=nrow(truth), ncol=Nsim)
HTse <- mod1se <- mod2se <- mod3se <- matrix(NA, nrow=nrow(truth), ncol=Nsim)

set.seed(1)
pums$id <- 1:nrow(pums)
keeps <- c()
for(sec in unique(pums$SECTOR)){
  for(st in unique(pums$FIPST)){
    tmp <- pums %>% filter(SECTOR==sec & FIPST==st)
    kid <- sample(tmp$id, 2)
    keeps <- c(keeps,kid)
  }
}

for(s in 1:Nsim){
  
  prob <- inclusionprobabilities(1/pums$TABWGT, ss)
  prob <- prob*ss/sum(prob)
  prob[keeps] <- 1
  Ind <- UPpoisson(prob)
  # Create sampled dataframe
  sample <- pums[as.logical(Ind),]
  sample$P <- prob[as.logical(Ind)]
  sample$W <- 1/ sample$P
  sample$scaledWGT <- (sample$W)*nrow(sample)/sum(sample$W)
  
  # Create direct estimate
  agg <- sample %>% group_by(SECTOR, FIPST) %>% summarize(HT_REC=horvitzThompson(y=RECEIPTS, pi=1/W)$pop_mean,
                                                        VAR_REC=horvitzThompson(y=RECEIPTS, pi=1/W, var_est=T, var_method = 'lin_HH')$pop_mean_var,
                                                        UW_REC=mean(RECEIPTS))
  
  HTests[,s] <- agg$HT_REC; UWests[,s] <- agg$UW_REC
  HTse[,s] <- sqrt(agg$VAR_REC)

  
  agg <- agg %>% left_join(predXT, by=c("SECTOR", "FIPST"))
  agg$VAR <- agg$VAR_REC/(agg$HT_REC+1)^2
  agg$VAR[agg$VAR==0] <- mean(agg$VAR)
  

  mod1 <- mod1Fit(Y=log(agg$HT_REC+1), X=cbind(1,log(agg$EMP),agg$ECOM), Var=agg$VAR, iter=1000, burn=500)
  mod1ests[,s] <- rowMeans(exp(mod1$Preds))
  mod1se[,s] <- apply(exp(mod1$Preds),1,sd)
  
  mod2 <- mod2Fit(Y=log(agg$HT_REC+1), 
                  X=cbind(1,log(agg$EMP),agg$ECOM), 
                  Psi1=model.matrix(~as.factor(agg$SECTOR)-1),
                  Psi2=model.matrix(~as.factor(agg$FIPST)-1),
                  Var=agg$VAR, iter=1000, burn=500)
  mod2ests[,s] <- rowMeans(exp(mod2$Preds))
  mod2se[,s] <- apply(exp(mod2$Preds),1,sd)
  
  sampX <- data.frame(SECTOR=sample$SECTOR,FIPST=sample$FIPST) %>% 
    left_join(agg%>%dplyr::select(SECTOR,FIPST,EMP,ECOM), by=c("SECTOR","FIPST"))
  predX <- cbind(1, log(pums$EMP+1), pums$ECOM)
  predPsi1 <- Matrix(model.matrix(~factor(pums$SECTOR, levels=levels(as.factor(sample$SECTOR)))-1))
  predPsi2 <- Matrix(model.matrix(~factor(pums$FIPST, levels=levels(as.factor(sample$FIPST)))-1))
  mod3 <- mod3Fit(Y=log(sample$RECEIPTS+1), 
                  X=cbind(1,log(sample$EMP+1),sample$ECOM), 
                  Psi1=model.matrix(~as.factor(sample$SECTOR)-1),
                  Psi2=model.matrix(~as.factor(sample$FIPST)-1),
                  wgt=sample$scaledWGT, 
                  predX=predX, 
                  predPsi1=predPsi1,
                  predPsi2=predPsi2,
                  predGrp=pums$Grp,
                  iter=1000, burn=500)
  mod3ests[,s] <- rowMeans(exp(mod3$Preds[,-1]))
  mod3se[,s] <- apply(exp(mod3$Preds[,-1]),1,sd)
  
  print(s)
  print(mean((HTests - truth$Mean)^2, na.rm=T)/10000000)
  print(mean((mod3ests - truth$Mean)^2, na.rm=T)/10000000)
}

simData <- list(truth, HTests, UWests, mod1ests, mod2ests, mod3ests, HTse, mod1se, mod2se, mod3se)
write_rds(simData, "sim_results_1_15_2020.rds")
## MSE
m1 <- mean((HTests - truth$Mean)^2, na.rm=T)
m5 <- mean((UWests - truth$Mean)^2, na.rm=T)
m2 <- mean((mod1ests - truth$Mean)^2, na.rm=T)
m3 <- mean((mod2ests - truth$Mean)^2, na.rm=T)
m4 <- mean((mod3ests - truth$Mean)^2, na.rm=T)
formatC(c(m1,m2,m3,m4,m5), format='e', digits=2)

## abs Bias
b1 <- mean(abs(rowMeans(HTests, na.rm=T) - truth$Mean), na.rm=T)
b5 <- mean(abs(rowMeans(UWests, na.rm=T) - truth$Mean), na.rm=T)
b2 <- mean(abs(rowMeans(mod1ests, na.rm=T) - truth$Mean), na.rm=T)
b3 <- mean(abs(rowMeans(mod2ests, na.rm=T) - truth$Mean), na.rm=T)
b4 <- mean(abs(rowMeans(mod3ests, na.rm=T) - truth$Mean), na.rm=T)
formatC(c(b1,b2,b3,b4,b5), format='e', digits=2)

plot(rowMeans(mod1se)/rowMeans(HTse))
plot(rowMeans(mod2se)/rowMeans(HTse))
plot(rowMeans(mod3se)/rowMeans(HTse), ylim=c(0,1))




