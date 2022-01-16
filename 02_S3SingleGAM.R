
# R script to replicate the single GAM analysis presented in S3 of
# Jacobson et al., Quantifying the response of Blainville's beaked whales 
# to US Naval sonar exercises in Hawaii
# Last updated 2022-01-16 by EKJ

m0 <- gam(DivePresent ~
            s(ID, bs="mrf", xt=list(nb=nb)) +
            s(Depth, bs="ts") +
            s(MaxRL, by=as.numeric(Sonar), bs="ts") +
            Phase +
            offset(logArea),
          family = binomial, data=m0.df, method="REML")

megagrid <- expand.grid(ID = unique(m0.df$ID), 
                        MaxRL=seq(100, 165, by= 1),
                        Phase = c("0", "A", "B")) 
megagrid$MaxRL[which(megagrid$Phase == "0" | megagrid$Phase == "A")] <- 0

megagrid <- left_join(m0.df[,c(21,20,15)], megagrid, by = "ID") %>%
  distinct() 

megagrid$Sonar <- 0
megagrid$Sonar[which(megagrid$Phase=="B")] <- 1

megagrid$m0.pred <- predict.gam(m0, newdata=megagrid, type = "response")

presonar <- filter(megagrid, Phase != "B") %>%
  pivot_wider(id_cols = c(ID), names_from = Phase, values_from = c(m0.pred))
withsonar <- filter(megagrid, Phase == "B")
m0.results <- left_join(withsonar, presonar, by = "ID")

m0.results$delta31 <- (m0.results$m0.pred-m0.results$"0")/m0.results$"0"
m0.results$delta32 <- (m0.results$m0.pred-m0.results$"A")/m0.results$"A"
m0.results$delta21 <- (m0.results$"A"-m0.results$"0")/m0.results$"0"

### Uncertainty Propagation

lp0 <- predict(m0, newdata=megagrid, type="lpmatrix")

m0.bs <- list()

for (i in 1:1000){
  
  nbeta0 <- rmvn(1, coef(m0), vcov(m0))
  megagrid$m0.pred <- as.numeric(megagrid$logArea + (lp0%*%nbeta0))
  megagrid$m0.exp <- m0$family$linkinv(megagrid$m0.pred)
  
  presonar <- filter(megagrid, Phase != "B") %>%
    pivot_wider(id_cols = c(ID), names_from = Phase, values_from = c(m0.pred, m0.exp))
  withsonar <- filter(megagrid, Phase == "B")
  
  m0.bs.results <- left_join(withsonar, presonar, by = "ID") 
  m0.bs.results$m0.pred_12 <- predict(m0, newdata = m0.bs.results)
  m0.bs.results$m0.exp_12 <- m0$family$linkinv(m0.bs.results$m0.pred_12)
  
  m0.bs.results$delta21 <- (m0.bs.results$m0.exp_A-m0.bs.results$m0.exp_0)/m0.bs.results$m0.exp_0
  m0.bs.results$delta31 <- (m0.bs.results$m0.exp_12-m0.bs.results$m0.exp_0)/m0.bs.results$m0.exp_0
  m0.bs.results$delta32 <- (m0.bs.results$m0.exp_12-m0.bs.results$m0.exp_A)/m0.bs.results$m0.exp_A
  
  m0.bs[[i]] <- m0.bs.results
  
  print(paste("Completed iteration", i))
}

m0.bootstrap <- as.data.frame(matrix(nrow = nrow(m0.results)*1000, ncol = 17))
names(m0.bootstrap) <- names(m0.bs[[1]])

for (i in 1:length(m0.bs)){
  
  start <- i*nrow(m0.results)-nrow(m0.results)+1
  end <- i*nrow(m0.results)
  m0.bootstrap[start:end, ] <- m0.bs[[i]]
}

m0.bootstrap.df <- m0.bootstrap

m0.d21 <- m0.bootstrap.df %>% group_by(MaxRL) %>%
  summarize("Q2.5" = quantile(delta21, probs=0.025),
            "Q50" = median(delta21),
            "Q97.5" = quantile(delta21, probs=0.975))

m0.d31 <- m0.bootstrap.df %>% group_by(MaxRL) %>%
  summarize("Q2.5" = quantile(delta31, probs=0.025),
            "Q50" = median(delta31),
            "Q97.5" = quantile(delta31, probs=0.975))

m0.d32 <- m0.bootstrap.df %>% group_by(MaxRL) %>%
  summarize("Q2.5" = quantile(delta32, probs=0.025),
            "Q50" = median(delta32),
            "Q97.5" = quantile(delta32, probs=0.975))

p.d31 <- ggplot(m0.d31) +
  geom_ribbon(aes(x=MaxRL, ymin=Q2.5, ymax=Q97.5), fill="grey") +
  geom_line(aes(x=MaxRL, y=Q50))+
  geom_rug(data=m0.df, aes(x=MaxRL))+
  ylab("Proportion change in P(GVP) rel. to Baseline") +
  xlab("MFAS Received Level (dB re. 1 \u00B5Pa)") +
  xlim(c(100, 165))+
  ylim(c(-1, 4))+
  theme_bw()

p.d32 <- ggplot(m0.d32) +
  geom_ribbon(aes(x=MaxRL, ymin=Q2.5, ymax=Q97.5), fill="grey") +
  geom_line(aes(x=MaxRL, y=Q50))+
  geom_rug(data=m0.df, aes(x=MaxRL))+
  ylab("Proportion change in P(GVP) rel. to Naval Training") +
  xlab("MFAS Received Level (dB re. 1 \u00B5 Pa)") +
  xlim(c(100, 165))+
  ylim(c(-1, 4))+
  theme_bw()

ggarrange(p.d32, p.d31)
