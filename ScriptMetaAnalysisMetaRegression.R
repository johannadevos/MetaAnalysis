#### Script written by Johanna de Vos, last updated on 14-02-2018.
#### It accompanies the following publication in Language Learning:
#### de Vos, Schriefers, Nivard & Lemhöfer (2018). A Meta-Analysis and Meta-Regression 
####   of Incidental Second Language Word Learning from Spoken Input.

## Read in library
library(metafor)

## If you don't have this library installed, follow this step:
#install.packages("metafor")

## Clear workspace
rm(list = ls())

## Set working directory
## NB: Refer to the directory on your computer where the text files with effect sizes and the covariance matrix are stored
setwd("C:/...")

## Read in data
## NB: read in the + or - 0.10 Effect sizes or Covariance matrix to conduct the sensitivity analysis
es <- read.table(file = "Effect sizes.txt", sep = "\t", header = TRUE, na.strings = "")
cov <- read.table(file = "Covariance matrix.txt", sep = "\t", header = FALSE, na.strings = "")

## Give each case its own ID and make it into a factor
es$ID <- 1:nrow(es)
es$ID <- as.factor(es$ID)

## Transform the covariance matrix from a dataframe into a matrix
cov <- as.matrix(cov)

## Extra: use peer-reviewed data only
#es_peer <- es[-which(es$peer_reviewed == "no"),]
#cov_peer <- cov[-which(es$peer_reviewed=="no"),-which(es$peer_reviewed=="no")]

## Run the null model (without any fixed or random effects)
mod0 <- rma.mv(Hedges_g, cov, data = es)
summary(mod0)

## Run Model 1a (without fixed effects, with a random intercept at the study level)
mod1a <- rma.mv(Hedges_g, cov, random = ~ 1 | study, data = es)
summary(mod1a)

## Run Model 1b (without fixed effects, with random intercepts at the study and sample level)
mod1b <- rma.mv(Hedges_g, cov, random = ~ 1 | study/ID, data = es)
summary(mod1b)

## Compare the models
anova(mod0, mod1a)
anova(mod1a, mod1b)

## Draw profile likelihood plots for Model 1b (see Appendix 4)
tiff("figures/Chapter 2 - Figure A.tiff", units="in", width=4, height=6, res=300)
par(mfrow=c(2,1))
profile(mod1b, sigma2=1)
profile(mod1b, sigma2=2)
dev.off()

## Calculate intra-class correlation (ICC) for Model 1b (the best-fitting Model 1)
round(mod1b$sigma2[1] / sum(mod1b$sigma2), 3)

## Default order of factor levels
es$age <- factor(es$age, levels = c("elementary school / kindergarten", "high school", "university"))
es$treatment <- factor(es$treatment, levels = c("audio", "audiovisual", "TBLT (-int)", "TBLT (+int)"))

## Run Model 2 (with five fixed effects, and with random intercepts at the study and sample level)
## Unreduced dataset
mod2 <- rma.mv(Hedges_g, cov, mods = ~ age + treatment + testing + gain_scores + control_group, random = ~ 1 | study/ID, data = es)
summary(mod2)

## Relevelling: change the order of factor levels
## Then rerun Model 2 to obtain significances for the other contrasts
es$age <- factor(es$age, levels = c("high school", "elementary school / kindergarten", "university"))
es$treatment <- factor(es$treatment, levels = c("audiovisual", "audio", "TBLT (-int)", "TBLT (+int)"))
es$treatment <- factor(es$treatment, levels = c("TBLT (-int)", "audio", "audiovisual", "TBLT (+int)"))

## Write the output of Model 2 to a file
output <- cbind(mod2$b, mod2$se, mod2$zval, mod2$pval, mod2$ci.lb, mod2$ci.ub)
colnames(output) <- cbind("Beta", "SE", "z", "p", "Lower bound", "Upper bound")
write.table(output, "Output Model 2 unreduced.txt", dec = ",")


#### Investigate influential cases

## Calculate Cook's distance
es$cookMod2 <- cooks.distance(mod2)

## Visualise outcomes
par(mfrow=c(1,1))
plot(es$cookMod2, ylab = "Cook's distance")

## Remove cases where Cook's distance > 1 (this only applies to Model 2)
largeCook2 <- es$cookMod2 > 0.85
es_cook2 <- es[-which(largeCook2 == TRUE),]
cov_cook2 <- cov[-which(largeCook2 == TRUE),-which(largeCook2 == TRUE)]

## In the sensitivity analysis, remove the same effect sizes as in the original analysis of the reduced dataset
#es_cook2 <- es[-(c(13, 103, 104)),]
#cov_cook2 <- cov[-(c(13, 103, 104)),-(c(13, 103, 104))]

## Rerun Model 2 with influential cases removed
mod2_cook <- rma.mv(Hedges_g, cov_cook2, mods = ~ age + treatment + testing + gain_scores + control_group, random = ~ 1 | study/ID, data = es_cook2)
summary(mod2_cook)

## Relevelling: change the order of factor levels
## Then rerun Model 2 to obtain significances for the other contrasts
es_cook2$age <- factor(es_cook2$age, levels = c("high school", "elementary school / kindergarten", "university"))
es_cook2$treatment <- factor(es_cook2$treatment, levels = c("audiovisual", "audio", "TBLT (-int)", "TBLT (+int)"))
es_cook2$treatment <- factor(es_cook2$treatment, levels = c("TBLT (-int)", "audio", "audiovisual", "TBLT (+int)"))
es_cook2$treatment <- factor(es_cook2$treatment, levels = c("TBLT (+int)", "audio", "audiovisual", "TBLT (-int)"))

## Go back to the default order of factor levels
es_cook2$age <- factor(es_cook2$age, levels = c("elementary school / kindergarten", "high school", "university"))
es_cook2$treatment <- factor(es_cook2$treatment, levels = c("audio", "audiovisual", "TBLT (-int)", "TBLT (+int)"))

## Draw profile likelihood plots for Model 2 (see Appendix 4)
tiff("figures/Chapter 2 - Figure B.tiff", units="in", width=4, height=6, res=300)
par(mfrow=c(2,1))
profile(mod2_cook, sigma2=1)
profile(mod2_cook, sigma2=2)
dev.off()

## Calculate intra-class correlation (ICC) for Model 2
round(mod2_cook$sigma2[1] / sum(mod2_cook$sigma2), 3)

## Write the output to a file
output_cook2 <- cbind(mod2_cook$b, mod2_cook$se, mod2_cook$zval, mod2_cook$pval, mod2_cook$ci.lb, mod2_cook$ci.ub)
colnames(output_cook2) <- cbind("Beta", "SE", "z", "p", "Lower bound", "Upper bound")
write.table(output_cook2, "Output Model 2 reduced.txt", dec = ",")


#### Extra (not in the article)

## Rerun Model 1 and Model 2 (unreduced) using maximum likelihood estimation (instead of restricted maximum likelihood)
## This is necessary for comparing models with a different number of fixed effects 
mod1b_ml <- rma.mv(Hedges_g, cov, random = ~ 1 | study/ID, data = es, method = "ML")
mod2_ml <- rma.mv(Hedges_g, cov, mods = ~ age + treatment + testing + gain_scores + control_group, random = ~ 1 | study/ID, data = es, method = "ML")

## Compare Model 1 to Model 2 (unreduced)
anova(mod1b_ml, mod2_ml)

## A comparison between Model 1 and Model 2 (reduced) is not possible, because they are computed on different size datasets