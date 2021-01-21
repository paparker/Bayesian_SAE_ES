library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mase)
library(ggthemes)
library(sampling)
library(sae)
source('R/saeFuncs.R')

## Read in direct estimates
directEsts <- read_rds('Data/DIRECT.rds')
directEsts$REC_SE[directEsts$REC_SE==0] <- max(directEsts$REC_SE)

## Fit Model 
set.seed(1)
mod2 <- mod2Fit(Y=log(directEsts$REC), 
                X=cbind(1,log(directEsts$EMP), directEsts$ECOM), 
                Psi1=model.matrix(~as.factor(directEsts$SECTOR)-1),
                Psi2=model.matrix(~as.factor(directEsts$FIPST)-1),
                Var=directEsts$REC_SE/(directEsts$REC^2), ## Using delta method for log transform
                iter=1000, burn=500)

## Generate output
modEst <- apply(exp(mod2$Preds),1,mean)
dirEst <- directEsts$REC
dirSE <- directEsts$REC_SE
modSE <- apply(exp(mod2$Preds),1,sd)
plotDF <- data.frame(dirEst, modEst, dirSE, modSE, Sector=directEsts$SECTOR, FIPST=directEsts$FIPST)

ggplot(plotDF)+
  facet_wrap(~Sector, scales='free')+
  geom_point(color='blue', alpha=0.5, aes(x=dirEst, y=modEst))+
  theme_classic()+
  geom_abline(slope=1, intercept=0, color='red')+
  xlab("Direct Estimate")+
  ylab("Model Estimate")+
  theme(axis.text.x = element_blank())
ggsave('pointEst.png')

ggplot(plotDF)+
  geom_violin(fill='blue', alpha=0.5, aes(x=Sector,y=modSE/dirSE))+
  facet_wrap(~Sector, scales='free')+
  theme_classic()+
  ylab("Ratio of Model to Direct Standard Errors")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave('seCompare.png')
