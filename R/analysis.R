library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mase)
library(ggthemes)
library(sampling)
library(sae)
source('R/saeFuncs.R')

setwd("/media/ryan/hdd/documents/Bayesian_SAE_ES/")
## Read in direct estimates
directEsts <- read_rds('Data/DIRECT.rds')
w <- which(directEsts$REC_SE==0)
directEsts <- directEsts[-w, ]
#directEsts$REC_SE[directEsts$REC_SE==0] <- max(directEsts$REC_SE)

## Fit Model
set.seed(314159)
Y = log(directEsts$REC)
X = cbind(1, log(directEsts$EMP), directEsts$ECOM)
sects <- sort(unique(directEsts$SECTOR))
# look at some linear model fits to see if the linear predictor is reasonable
idx <- 20
tmpFit <- lm(Y ~ X - 1, subset = (directEsts$SECTOR == sects[idx]))
summary(tmpFit)
Y_hat <- tmpFit$fitted.values
plot(Y_hat, Y[directEsts$SECTOR == sects[idx]], main = sects[idx])
abline(0, 1, col = 'red')
plot(exp(Y_hat), exp(Y[directEsts$SECTOR == sects[idx]]), main = sects[idx])
abline(0, 1, col = 'red')

mod2 <- mod2Fit(Y = log(directEsts$REC),
                X = cbind(1, log(directEsts$EMP), directEsts$ECOM),
                Psi1 = model.matrix(~as.factor(directEsts$SECTOR) - 1),
                Psi2 = model.matrix(~as.factor(directEsts$FIPST) - 1),
                # there was a problem here
                Var = (directEsts$REC_SE^2/directEsts$REC^2),          ## Using delta method for log transform
                iter = 10000, burn = 5000)

## Generate output
modEst <- apply(exp(mod2$Preds),1,mean)
dirEst <- directEsts$REC
plot(modEst, dirEst)
abline(0, 1, col = 'red')
# this looks good!!

dirSE <- directEsts$REC_SE
modSE <- apply(exp(mod2$Preds), 1, sd)
plot(dirSE, modSE)
abline(0, 1, col = 'red')

plotDF <- data.frame(dirEst, modEst, dirSE, modSE, Sector = directEsts$SECTOR, FIPST = directEsts$FIPST)

ggplot(plotDF)+
  facet_wrap(~Sector, scales='free')+
  geom_point(color='blue', alpha=0.5, aes(x=dirEst, y=modEst))+
  theme_classic()+
  geom_abline(slope=1, intercept=0, color='red')+
  xlab("Direct Estimate")+
  ylab("Model Estimate")+
  theme(axis.text.x = element_blank())
ggsave('pointEst.png')
# this looks better; still a few outliers, but overall much better

ggplot(plotDF)+
  geom_violin(fill='blue', alpha=0.5, aes(x=Sector,y=modSE/dirSE))+
  facet_wrap(~Sector, scales='free')+
  theme_classic()+
  ylab("Ratio of Model to Direct Standard Errors")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave('seCompare.png')
# this looks much more reasonable; a few outliers in sector 51, 99, 
#  and something strange happening in sector 22


idx <- 6
plot(exp(mod2$Preds[idx, ]), main = idx,
     ylim = c(min(exp(Y[idx]), min(exp(mod2$Preds[idx, ]))), max(exp(Y[idx]), max(exp(mod2$Preds[idx])))))
abline(exp(Y[idx]), 0, col = 'red')
dirSE[idx]
# 98890.75
modSE[idx]
# 67.64231
modSE[idx]/dirSE[idx]
# 0.0006840104
# this doesn't make sense

# let's try using the sae package industry by industry
idx1 <- 1
dirEst <- directEsts[directEsts$SECTOR == sects[idx1], ] %>%
  dplyr::mutate(LOG_EMP = log(EMP)) %>%
  dplyr::mutate(LOG_REC = log(REC)) %>%
  dplyr::mutate(LOG_PAY = log(PAY)) %>%
  dplyr::mutate(LOG_REC_SE = REC_SE/REC) %>%
  dplyr::mutate(LOG_REC_VAR = (REC_SE/REC)^2)

saeFH <- sae::eblupFH(LOG_REC ~ LOG_EMP + LOG_PAY, vardir = LOG_REC_VAR, data = as.data.frame(dirEst))
Y_hat <- saeFH$eblup
plot(dirEst$LOG_REC, Y_hat)
abline(0, 1, col = 'red')
saeFH_MSE <- sae::mseFH(LOG_REC ~ LOG_EMP + LOG_PAY, vardir = LOG_REC_VAR, data = as.data.frame(dirEst))
Y_hat_MSE <- saeFH_MSE$mse
plot(dirEst$LOG_REC_VAR, Y_hat_MSE)
abline(0, 1, col = 'red')
summary(Y_hat_MSE/dirEst$LOG_REC_VAR)
