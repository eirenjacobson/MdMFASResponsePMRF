# R script to replicate the multi-GAM analyses presented in 
# Jacobson et al., Quantifying the response of Blainville's beaked whales 
# to US Naval sonar exercises in Hawaii
# Last updated 2021-11-03 by EKJ

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(scam)
library(ggpubr)

# Load the PMRF dataset
load("./Data/PMRF_AllSCCs.RData")

# Load the list of all HPs with their neighbors
# this is from the tessellation made using deldir
load("./Data/PMRF_HPNeighbors.RData")

# Convert the info in neighbors into the indexed format needed for the MRF
nb <- HPNeighbors
for (i in 1:length(nb)){
  nb[[i]] <- which(names(HPNeighbors) %in% HPNeighbors[[i]])
}

# baseline data where no training or mfas was present
baseline <- filter(PMRF_SCCData, Period=="00") 
# training where only training, no mfas was present
training <- filter(PMRF_SCCData, Period=="02")
# sonar when both training and mfas were present
sonar <- filter(PMRF_SCCData, Period == "12")
# within sonar dataset, if sonar RL was estimated to be 0 mark as no sonar
sonar$Sonar[which(sonar$Sonar==1 & sonar$MaxRL == 0)] <- 0
# filter out times when no sonar was present (e.g., in between exercises)
sonar <- filter(sonar, MaxRL>0)

# bind together the baseline, training, and sonar datasets for use when
# fitting a single GAM (M0)
m0.df <- rbind.data.frame(baseline, training, sonar)

## M1: Modelling the pre-activity probability of dive detection

# First, we model the spatial distribution of beaked whale detections 
# during the baseline period (when no naval training activity nor sonar 
# were present) with a Markov Random Field.  We include a covariate for 
# depth and use log(Area) of the tessellation tile around each hydrophone 
# as an offset. 

m1 <- gam(DivePresent ~
            s(ID, bs="mrf", xt=list(nb=nb)) +
            s(Depth, bs="ts") +
            offset(logArea),
          family = binomial, data=baseline,
          method="REML")

m1.predgrid <- baseline %>%
  group_by(ID, Depth, logArea) %>%
  summarize(PDive = sum(DivePresent)/length(DivePresent))

m1.predgrid$m1.pred <- predict.gam(m1, newdata=m1.predgrid, type="response")

ggplot(m1.predgrid) +
  geom_point(aes(x=PDive, y = m1.pred)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Observed probability of a GVP detection") +
  ylab("Predicted probability of a GVP detection") +
  theme_bw()

# Then, we use this model to predict the mean values expected by the 
# baseline model (on the logit scale) onto the data from when naval 
# training was present to be used as an offset 

# pull out the unique combinations of predictors for m1
m1.predgrid <- distinct(baseline, ID, .keep_all=TRUE) %>% 
  select(ID, Depth, logArea)

# predict expected values for each combo
m1.predgrid$m1.pred <- as.numeric(predict.gam(m1, newdata=m1.predgrid))

# stick these predicted values onto the training dataset
training <- left_join(training, m1.predgrid[,c(1,4)], by = c("ID"))
training$m1.pred <- as.numeric(training$m1.pred)

## M2: Modelling the effect of Naval training activity

# We then fit a model to the data from when naval training was present, 
# using the expected value from the baseline model (M1) as an offset 

m2 <- gam(DivePresent ~ 
            offset(m1.pred), 
          family = binomial, data = training, method = "REML")

# As with M1 to M2, we now use M2 to predict the expected values 
# at times when sonar was present 

# predict expected P(GVP) using M1 and M2
sonar$m1.pred <- as.numeric(predict.gam(m1, newdata=sonar))
sonar$m2.pred <- as.numeric(predict.gam(m2, newdata=sonar))

## M3: Modelling the effect of hull-mounted MFAS

m3 <- scam(DivePresent ~ s(MaxRL, bs="mpd") + offset(m2.pred) - 1,
           family=binomial, data=sonar)

predgrid <- distinct(baseline, ID, .keep_all=TRUE) %>% 
  select(ID, Depth, logArea)

megagrid <- expand.grid(ID = predgrid$ID, 
                        MaxRL=seq(100, 165, by= 1))

megagrid <- left_join(predgrid, megagrid, by = "ID")
megagrid$m1.pred <- as.numeric(predict(m1, newdata=megagrid))
megagrid$m2.pred <- predict(m2, newdata=megagrid)
megagrid$m3.pred <- predict(m3, newdata=megagrid)

megagrid$m1.exp <- m1$family$linkinv(megagrid$m1.pred)
megagrid$m2.exp <- m2$family$linkinv(megagrid$m2.pred)
megagrid$m3.exp <- m3$family$linkinv(megagrid$m3.pred)

megagrid$delta1 <- (megagrid$m2.exp-megagrid$m1.exp)/megagrid$m1.exp # ships rel. to baseline
megagrid$delta2 <- (megagrid$m3.exp-megagrid$m2.exp)/megagrid$m2.exp # sonar rel. to ships
megagrid$delta3 <- (megagrid$m3.exp-megagrid$m1.exp)/megagrid$m1.exp # sonar rel. to baseline

megasummary <- megagrid %>%
  group_by(MaxRL) %>%
  summarize(delta3 = mean(delta3),
            delta2 = mean(delta2))

ggplot(megasummary, aes(x=MaxRL, y=delta3))+
  geom_line()+
  ylab("Change Relative to Baseline") +
  ylim(c(-1, 0))+
  xlim(100, 165)+
  theme_bw()

ggplot(megasummary, aes(x=MaxRL, y=delta2))+
  geom_line()+
  ylab("Change Relative to Training") +
  ylim(c(-1, 0))+
  xlim(100, 165)+
  theme_bw()

## Uncertainty Propagation

m1.predgrid.bs <- distinct(baseline, ID, .keep_all=TRUE) %>% 
  select(ID, Depth, logArea)
lp1.pg <- predict.gam(m1, newdata=m1.predgrid.bs, type = "lpmatrix")

# get the linear predictors for m2
lp2.s <- predict.gam(m2, newdata=sonar, type="lpmatrix")

# set up prediction data frames
megagrid.bs <- expand.grid(ID = predgrid$ID, 
                        MaxRL=seq(100, 165, by= 1))
megagrid.bs <- left_join(m1.predgrid.bs, megagrid.bs, by = "ID")

# fill in preds for m2 and m3 so we can get the matrix of linear predictors
megagrid.bs$m1.pred <- as.numeric(predict.gam(m1, newdata=megagrid.bs))
megagrid.bs$m2.pred <- as.numeric(predict.gam(m2, newdata=megagrid.bs))

# pull out the lp matrices corresponding to the prediction grid
lp1m <- predict(m1, newdata=megagrid.bs, type="lpmatrix")
lp2m <- predict(m2, newdata=megagrid.bs, type="lpmatrix")
lp3m <- predict(m3, newdata=megagrid.bs, type="lpmatrix")

##bootstrap.df <- data.frame()
bootstrap.list <- list()

set.seed(20211025)

for (i in 1:5001){
  
  if(i %in% c(1001, 2001, 3001, 4001, 5001)){
    bootstrap.mat <- do.call(rbind, bootstrap.list)
    save(bootstrap.mat, file = paste0("./Data/bootstrapM1M2M3_", i, ".RData"))
    bootstrap.list <- list()
  }
  
  # refilter the data to make sure no prev offsets are accidentally used
  training.bs <- training[,1:7]
  sonar.bs <- sonar[,1:7]
  
  # pull a random draw of the coefficients of m1
  nbeta1 <- rmvn(1, coef(m1), vcov(m1))
  
  # fill the m1.predgrid slot with the random draw on logit scale
  m1.predgrid.bs$m1.pred <- as.numeric(m1.predgrid.bs$logArea+lp1.pg%*%nbeta1)
  
  # now stick that random draw onto the other datasets to be used as an offset later
  training.bs <- left_join(training.bs, m1.predgrid.bs[,c(1,4)], by = c("ID"))
  training.bs$m1.pred <- as.numeric(training.bs$m1.pred)
  sonar.bs <- left_join(sonar.bs, m1.predgrid.bs[,c(1,4)], by = c("ID"))
  sonar.bs$m1.pred <- as.numeric(sonar.bs$m1.pred)
  
  # Model 2: Training data with random draw from predicted baseline values as offset
  
  # Fit M2 using random draw from baseline as offset  
  m2.bs <- gam(DivePresent ~ 
              offset(m1.pred), 
            family = binomial, data = training, method = "REML")
  
  # Model 3: Sonar data with predicted ship values as offset
  
  # pull a random draw of the coefficients of m2
  nbeta2 <- rmvn(1, coef(m2.bs), vcov(m2.bs))
  
  # fill the m2.pred slot with the random draw on logit scale
  sonar.bs$m2.pred <- as.numeric(sonar$m1.pred+lp2.s%*%nbeta2)
  
  # fit M3
  m3.bs <- scam(DivePresent ~ s(MaxRL, bs="mpd", m=2) + offset(m2.pred) - 1,
             family=binomial, data=sonar)
  
  nbeta3 <- rmvn(1, m3.bs$coefficients.t, m3.bs$Vp.t)

  # make predictions using new models 
  megagrid.bs$m1.pred <- as.numeric(megagrid.bs$logArea+lp1m%*%nbeta1)
  megagrid.bs$m2.pred <- as.numeric(megagrid.bs$m1.pred+lp2m%*%nbeta2)
  megagrid.bs$m3.pred <- as.numeric(megagrid.bs$m2.pred + lp3m%*%nbeta3)
  
  megagrid.bs$m1.exp <- m1$family$linkinv(megagrid.bs$m1.pred)
  megagrid.bs$m2.exp <- m2$family$linkinv(megagrid.bs$m2.pred)
  megagrid.bs$m3.exp <- m3$family$linkinv(megagrid.bs$m3.pred)
  
  # calculate change
  megagrid.bs$delta1 <- (megagrid.bs$m2.exp-megagrid.bs$m1.exp)/megagrid.bs$m1.exp # training rel. to baseline
  megagrid.bs$delta2 <- (megagrid.bs$m3.exp-megagrid.bs$m2.exp)/megagrid.bs$m2.exp # sonar rel. to training
  megagrid.bs$delta3 <- (megagrid.bs$m3.exp-megagrid.bs$m1.exp)/megagrid.bs$m1.exp # sonar rel. to baseline
  
  megagrid.bs$Iteration <- i
  
  # attach this iteration to the main df
  bootstrap.list[[i]] <- as.matrix(select(megagrid.bs, c(Iteration, MaxRL, delta1, delta2, delta3)))
  
  print(paste("Completed interation", i))
  
}

load("./Data/bootstrapM1M2M3_1001.RData")
bs1 <- bootstrap.mat
load("./Data/bootstrapM1M2M3_2001.RData")
bs2 <- bootstrap.mat
load("./Data/bootstrapM1M2M3_3001.RData")
bs3 <- bootstrap.mat
load("./Data/bootstrapM1M2M3_4001.RData")
bs4 <- bootstrap.mat
load("./Data/bootstrapM1M2M3_5001.RData")
bs5 <- bootstrap.mat

bootstrap.df <- as.data.frame(rbind(bs1, bs2, bs3, bs4, bs5))
rm(bs1, bs2, bs3, bs4, bs5)

# training relative to baseline
d21 <- bootstrap.df %>%
  summarize("Q2.5" = quantile(delta1, probs=0.025),
            "Q50" = median(delta1),
            "Q97.5" = quantile(delta1, probs=0.975))

# sonar relative to training
d32 <- bootstrap.df %>%
  group_by(MaxRL) %>%
  summarize("Q2.5" = quantile(delta2, probs=0.025),
            "Q50" = median(delta2),
            "mean" = mean(delta2),
            "Q97.5" = quantile(delta2, probs=0.975))

# sonar relative to baseline
d31 <- bootstrap.df %>%
  group_by(MaxRL) %>%
  summarize("Q2.5" = quantile(delta3, probs=0.025),
            "Q50" = median(delta3),
            "mean" = mean(delta3),
            "Q97.5" = quantile(delta3, probs=0.975))

p31 <- ggplot(d31) +
  geom_ribbon(aes(x=MaxRL, ymin=Q2.5, ymax=Q97.5), fill="grey") +
  geom_line(aes(x=MaxRL, y=Q50))+
  geom_rug(data=sonar, aes(x=MaxRL))+
  ylab("Proportion change in P(GVP) rel. to Baseline") +
  xlab("MFAS Received Level (dB re 1 \u00B5Pa)") +
  xlim(c(100, 165))+
  ylim(c(-1, 0.1))+
  theme_bw()

p32 <- ggplot(d32) +
  geom_ribbon(aes(x=MaxRL, ymin=Q2.5, ymax=Q97.5), fill="grey") +
  geom_line(aes(x=MaxRL, y=Q50))+
  geom_rug(data=sonar, aes(x=MaxRL))+
  ylab("Proportion change in P(GVP) rel. to Naval Activity") +
  xlab("MFAS Received Level (dB re 1 \u00B5Pa)") +
  xlim(c(100, 165))+
  ylim(c(-1, 0.1))+
  theme_bw()

ggarrange(p31, p32)
