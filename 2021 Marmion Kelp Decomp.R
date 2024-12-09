#################################################################
# 2021 MARMION IN SITU KELP DECOMPOSITION   EXPERIMENT
# Taylor Simpkins: taylor.simpkins@research.uwa.edu.au ; taylorsimpkinsAUS@gmail.com
# UNIVERSITY OF WESTERN AUSTRALIA - Wernberg Lab
# Co-authors: Karen Filbee-Dexter, Mirjam Van Der Mheen, Albert Pessarrodona, Thomas Wernberg
#################################################################
citation ()

summary(decomp_2021)

setwd("~/Desktop/R_litterbag")
getwd()

attach (decomp_2021)
library(nlme)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(car)
library(MASS)
library(patchwork)
library(magrittr)
library(purrr)
library(cowplot) #theme_half_open
library(grDevices) #times new roman
library(simplecolors) #ggplot extra colors


decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$depth <- factor(decomp_2021$depth)
decomp_2021$species <- factor(decomp_2021$species)
decomp_2021$site <- factor(decomp_2021$site)
decomp_2021$gC <- ((decomp_2021$C/100) * decomp_2021$dw)
decomp_2021$gC
decomp_2021$D <- (((1.01-decomp_2021$rb)/decomp_2021$time)*100)
# decomp_2021$D <-  subset(decomp_2021, time == "14")
decomp_2021$D 
#decomp_2021$time <- factor(decomp_2021$time)
head(decomp_2021)

# visualise gC (%C x remaining biomass(g))
decomp_2021 %>% 
  ggplot(aes(time, gC, colour = species)) +
  geom_point() +
  facet_wrap(~species) +
  theme_bw()
# visualise remaining biomass (proportion of initial)
decomp_2021 %>% 
  ggplot(aes(time, rb, colour = species)) +
  geom_point() +
  facet_wrap(~species) +
  theme_bw()

decay_decomp_21 <- subset(decomp_2021, subset = time > 1)
decay_decomp_21$time <- factor(decay_decomp_21$time)
decay_decomp_21
decay_decomp_21 %>% 
  ggplot(aes(time, D, fill = species)) +
  geom_boxplot()+
  theme_bw()

# check for correlation between decay rate and %C x depth for each species
decay_decomp_21 %>% 
  filter(time == "14") %>% 
  ggplot(aes(C, D, group = depth, colour = depth)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~species) +
  theme_bw()

DCmod <- lm(C~D*depth, data=decay_decomp_21, species=="Ecklonia radiata")
summary(DCmod)
DCmod <- lm(C~D*depth, data=decay_decomp_21, species=="Scytothalia dorycarpa")
summary(DCmod)
# no effect of C on D or depth on regressions for either species

# check for correlation between decay rate and CN x depth for each species
decay_decomp_21 %>% 
  filter(time == "14") %>% 
  ggplot(aes(C.N, D, group = depth, colour = depth)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~species) +
  theme_bw()
# no effect of C on D or depth on regressions for either species

# visualise potential correlation between C:N and D
decay_decomp_21 %>% 
  ggplot(aes(C.N, D, colour=species)) +
  geom_point() +
  theme_bw()


#################################################################
########## PREDICTIVE MODEL FOR NEG EXP DECAY ###################
#################################################################

# H1: Kelp decomp follows a negative exponential decay model -> gnls
# H2 (NULL): No significant difference in percent remaining biomass (rb) between depths - summary of gnls model ()
# H3 (NULL): No significant difference in decay rates of species (m4)
# output table: k + A values

str(decomp_2021)
# SSasymp is a self starting function
# C is the area above the asymptote
# A is the area below the asymptote
# I (y int) should be close to 1 (100%)


###################################################################
# MODEL1 - %rb no factors ###### INTERSPECIFIC KELP DECOMP? #######

# m1 to determine start value and A for m2
m1 <- nls(rb ~ SSasymp(time, A, I, logk), data = decomp_2021) # m1 = all data: 10, 20, and 50m depths for both species
summary(m1)
# Formula: rb ~ SSasymp(time, A, I, logk)
# Parameters:
#       Estimate Std. Error t value Pr(>|t|)    
# A     0.11318    0.03364   3.364 0.000931 ***
# I     0.99131    0.05659  17.518  < 2e-16 ***
# logk -2.34805    0.12503 -18.780  < 2e-16 ***


# MODEL2 - %rb with species as factor ##############################

# NEW MODEL with values from m1 for start for each species
# write out the entire eqn: simplify by removing I since it is 1
m2 <- nls(rb ~ exp(k[species] * time),
          start = list(k = c(-exp(-2.33516),-exp(-2.33516))),
          data = decomp_2021)
summary(m2)
# Formula: rb ~ exp(k[species] * time)
# Parameters:
#       Estimate Std. Error t value Pr(>|t|)    
# k1 -0.058245   0.003844  -15.15   <2e-16 *** (ECKLONIA)
# k2 -0.089354   0.006454  -13.84   <2e-16 *** (SCYTOTHALIA)
plot(resid(m2)~decomp_2021$time)
abline(0,0)
hist(resid(m2))
confint(m2)


# MODEL3 - include A (residual carbon?) ############################

m3 <- nls(rb ~ exp(k[species] * time) + A[species],
          start = list(k = c(-exp(-2.33516),-exp(-2.33516)),
                       A = c(0.11318,0.11318)),
          data = decomp_2021)
summary(m3)
# Parameters:
#     Estimate Std. Error t value Pr(>|t|)    
# k1 -0.07602    0.01334  -5.698 4.68e-08 ***
# k2 -0.11419    0.01677  -6.808 1.32e-10 ***
# A1  0.08185    0.05271   1.553   0.1222    
# A2  0.07138    0.03597   1.985   0.0487 * 
plot(resid(m3)~decomp_2021$time)
abline(0,0)
hist(resid(m3))
# heteroscedastic


# MODEL4 - weights to account for interspecific variation in rb and lower n for final timepoint
decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata", "Scytothalia dorycarpa"))
m4 <- gnls(rb ~ exp(k * time) + A,
           start = list(k = c(-exp(-2.33516),-exp(-2.33516)),
                        A = c(0.11318,0.11318)),
           params = list(k ~ species,
                        A ~ species),
           weights = varExp(form = ~time|species),
           data = decomp_2021)
plot(resid(m4, type = "normalized") ~ time, data = decomp_2021)
abline(0,0) # PASS visual homogeneity test: weighting for time and species adjusts for heteroscedasticity at the final timepoint
hist(resid(m4, type = "normalized")) # PASS visual normality test 

# Assumptions of parametric tests met >> Test for interspecific variation in rb over time (decomposition rate)
# summary.gls function generates a tTable of summary statistics:
#   a matrix with columns Value, Std. Error, t-value, and p-value 
#   representing respectively the coefficients estimates, their approximate standard errors, the ratios between the estimates and their standard errors, 
#   and the associated P-VALUE under a t approximation. Rows correspond to the different coefficients.
# !!! USE SUMMARY ONLY FOR VALUE ESTIMATES, pvals are based on T-tests... NEED ANOVA
summary(m4) #ECKLONIA
intervals(m4)
# Coefficients:
#                                       Value  Std.Error   t-value p-value
# k.(Intercept)                  -0.07534915 0.01546449 -4.872397  0.0000 (kval for ECKLONIA ***)
# k.speciesScytothalia dorycarpa -0.04177590 0.01862877 -2.242548  0.0261 (kval for ECKLONIA slower than SCYTO *)
# A.(Intercept)                   0.07988614 0.06106027  1.308316  0.1924 (A not significant > no residual tissue (C))
# A.speciesScytothalia dorycarpa  0.00458024 0.06282662  0.072903  0.9420 (A for SCYTO not different from ECKLONIA)
# REVERSE ORDER OF SPECIES TO GET K VAL FOR SCYTO
decomp_2021$species <- factor(decomp_2021$species, levels = c("Scytothalia dorycarpa",
                                                              "Ecklonia radiata"))
m4 <- gnls(rb ~ exp(k * time) + A,
           start = list(k = c(-exp(-2.33516),-exp(-2.33516)),
                        A = c(0.11318,0.11318)),
           params = list(k ~ species,
                         A ~ species),
           weights = varExp(form = ~time|species),
           data = decomp_2021)
summary(m4) # SCYTOTHALIA
intervals(m4)
# Coefficients:
#                                 Value  Std.Error    t-value p-value
# k.(Intercept)             -0.11712637 0.01038672 -11.276549  0.0000 (kval for SCYTOTHALIA ***)
# k.speciesEcklonia radiata  0.04178092 0.01862868   2.242828  0.0261
# A.(Intercept)              0.08446773 0.01479273   5.710084  0.0000 (A significant *** > residual tissue (C))
# A.speciesEcklonia radiata -0.00459537 0.06282941  -0.073140  0.9418 (Ecklonia A borderline *)

# ANOVA
ggplot(decomp_2021, aes(x = (rb))) +   
  geom_histogram() + 
  facet_wrap(~ depth*species, ncol = 1) 
library(lme4)
fm3<-lmer(rb ~ depth*species*time+(1|site), data=decomp_2021)
#this is if you need to better fit your data's distribution. 
fm3<-glmer(rb ~ depth+species+Time+(1|Tank), family=Gamma(link = "identity" ),data=decomp_2021)
           
           summary(fm3)
           anova(fm3)
          
           
           #look at model fit
           plot(fitted(fm3), residuals(fm3), xlab = "Fitted Values", ylab = "Residuals")
           qqnorm(residuals(fm3));qqline(residuals(fm3))
           #ranova(fm3)

decomp_2021$time <- factor(decomp_2021$time)
model <- aov(rb~species*depth*time, data=decomp_2021)
summary(model)
#                 Df Sum Sq Mean Sq F value  Pr(>F)   
# species         1  0.104  0.1037   1.336 0.24930   
# depth           2  0.900  0.4500   5.797 0.00362 **
# species:depth   2  0.637  0.3184   4.102 0.01808 * 
# Residuals     184 14.283  0.0776 
TukeyHSD(model)
# $depth
# diff         lwr         upr     p adj
# 20-10 -0.1498887 -0.25778569 -0.04199165 0.0035123
# 50-10 -0.0420930 -0.17963525  0.09544926 0.7501227
# 50-20  0.1077957 -0.02182372  0.23741506 0.1238587
# Scytothalia dorycarpa:20-Scytothalia dorycarpa:10 -0.22112532 -0.40260344 -0.03964720 0.0073577
# Ecklonia radiata:50-Scytothalia dorycarpa:20       0.26347826  0.04039218  0.48656434 0.0105045
# Scytothalia dorycarpa:50:14-Ecklonia radiata:50:14      -0.007870636 0.0337862
# Scytothalia dorycarpa:20:14-Scytothalia dorycarpa:10:14  0.349351587 0.9853997
# Ecklonia radiata:50:14-Scytothalia dorycarpa:10:14       0.482129364 0.0319189

# Scytothalia dorycarpa:10:14-Ecklonia radiata:10:14      -0.079502924
# Scytothalia dorycarpa:20:14-Ecklonia radiata:20:14      -0.104027778
# Scytothalia dorycarpa:50:14-Ecklonia radiata:50:14      -0.244444444
# Scytothalia dorycarpa:10:31-Ecklonia radiata:10:31      -0.090000000
# Scytothalia dorycarpa:20:31-Ecklonia radiata:20:31       0.106563572

model <- aov(p~species*depth*time, data=decomp_2021)
summary(model)

model <- aov(rb~species*time, data=decomp_2021)
summary(model)


decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata", "Scytothalia dorycarpa"))
# Use Log Likelihood to show goodness of fit of model 4 to other models (Higher value = better model)
logLik(m1) #33.51
logLik(m2) #37.80
logLik(m3) #39.97
logLik(m4) #65.77

AIC(m1)
AIC(m2)
AIC(m3)
AIC(m4)

# # Calculate and compare AIC values
# install.packages("AICcmodavg")
# library(AICcmodavg)
# models <- list(m1e, m2e, m3e)
# mod.names <- c('species', 'depth')
# aictab(cand.set = models)#, modnames=mod.names)

# FIGURE 2.A ################################################
decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata", "Scytothalia dorycarpa"))
# Start by calculating means and SE for each species rb
decomp_summary_21 <- ddply (decomp_2021, c("species", "time"), summarise,
                            N    = sum(!is.na(rb)),
                            mean = mean(rb, na.rm=TRUE),
                            median = median(rb, na.rm=TRUE),
                            sd   = sd(rb, na.rm = TRUE), 
                            se   = sd / sqrt(N))
decomp_summary_21
# make a fit line for decay estimated by model 4
decomp_new <- data.frame(time = rep(seq(0, 50, by = 1), 2),
                         species = c(rep("Ecklonia radiata", 51), rep("Scytothalia dorycarpa", 51)))
decomp_new$species <- factor(decomp_new$species)
decomp_new$fit <- predict(m4, newdata =  decomp_new)
decomp_2021$time <- as.factor(decomp_2021$time)
#decomp_new$time <- numeric(decomp_new$time)
decomp_new
# decomp_new %<>% mutate (low = fit * 0.75, high = fit * 1.25, index = seq(1,102))
# CONFIDENCE INTERVALS WORKFLOW: https://stackoverflow.com/questions/32459480/r-confidence-bands-for-exponential-model-nls-in-basic-graphics
# ci = outer(decomp_new$fit, c(outer(decomp_summary_21$se, c(-1,1), '*'))*1.96, '+')
# ii = order(decomp_new$time)
# geom_(decomp_new[ii,], plot(time, fit, ylim=range(ci), type='l'))
# matlines(decomp_new[ii,'time'], ci[ii,], lty=2, col=1)

# time must be continuous
decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$species <- factor(decomp_2021$species)

particle_MAR$time <- numeric(particle_MAR$time)
# Plot rb means +-se for each species with decomp_new$fit
ggplot() + geom_point(data = decomp_2021, aes(time, rb, color=species), na.rm=TRUE, position=position_dodge(width=2), size = 1.5, alpha = 0.2) +
  # geom_smooth(data = decomp_2021, aes(time, rb, color = species)) +
  geom_line(data = decomp_new, aes(time, fit, color=species), size=1, alpha=0.6) +
  # geom_ribbon(data = decomp_new, aes(time, fit, shading=species, ymin = low, ymax = high), alpha = 0.1) +
  # geom_point(data = particle_MAR, aes(time, export,color = depth), size = 2.5, alpha = 0.5) +
  geom_hline(yintercept=c(0.0), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.25), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.75), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(1.0), linetype="dotted", alpha = 0.5) +
  # geom_vline(xintercept=c(11.3), linetype="dotted", alpha = 0.5) +
  # geom_vline(xintercept=c(7.6), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = decomp_summary_21, aes(time, mean, ymin = mean - (se), ymax = mean + (se)), na.rm = TRUE, width = 0.1) +
  geom_point(data = decomp_summary_21, aes(time, mean, fill=species),colour="black", na.rm = TRUE, pch = 21, size = 3) + # means added from summary stats
  #scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  ylim(0, 1) +
  labs(x="Days (in situ)", y="Proportion of initial carbon (dw)",  title="All species and depths") +  #title= "Interspecific kelp decomposition") +
  theme_half_open()+
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
ggsave("MAR_decomp_plot.pdf", width = 4, height = 4, units = "in")
  



###################################################################
# MODEL1E - ecklonia %rb with depth ################################

decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$depth <- factor(decomp_2021$depth)
decomp_2021$species <- factor(decomp_2021$species)
decomp_2021$site <- factor(decomp_2021$site)
decomp_2021$gC <- (decomp_2021$C * decomp_2021$dw)/100
decomp_2021$D <- ((1-decomp_2021$rb)/decomp_2021$time)
ecklonia_decomp_21 <- subset(decomp_2021, species == "Ecklonia radiata")
#calculate mean initial gC for E.radiata
mean(1.4028406, 1.4260579, 1.6029657) #[1] 1.402841
ecklonia_decomp_21$prC <- (ecklonia_decomp_21$gC/1.402841)
ecklonia_decomp_21$prC[is.na(ecklonia_decomp_21$prC)] <- 0
ecklonia_decomp_21$prC #column created for percent remaining carbon

#ecklonia_decomp_21$time <- numeric(ecklonia_decomp_21$time)

m1e <- nls(rb ~ SSasymp(time, A, I, logk), 
           data = ecklonia_decomp_21)
summary(m1e)
# Formula: rb ~ SSasymp(time, A, I, logk)
# Parameters:
#       Estimate Std. Error t value Pr(>|t|)    
# A     0.09029    0.08073   1.119    0.266    
# I     0.94357    0.14669   6.433 6.12e-09 ***
# logk -2.69769    0.31144  -8.662 1.88e-13 ***

# NEW MODEL with values from m1 for start for each depth
# write out the entire eqn: simplify by removing I since it is 1
m2e <- nls(rb ~ exp(k[depth] * time),
           start = list(k = c(-exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
           data = ecklonia_decomp_21)
summary(m2e)
# Formula: biomass ~ exp(k[depth] * time)
# Parameters:
# k1 -0.080187   0.012312  -6.513 4.26e-09 ***
# k2 -0.049351   0.006023  -8.194 1.74e-12 ***
# k3 -0.053174   0.008980  -5.921 5.89e-08 ***

m3e <- nls(rb ~ exp(k[depth] * time) + A[depth],
           start = list(I = c(0.94357,0.94357,0.94357),
                        k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769)),
                        A = c(0.09029,0.09029,0.09029)),
           data = ecklonia_decomp_21)
summary(m3e)
# Formula: biomass ~ exp(k[depth] * time) + A[depth]
# Parameters:
#       Estimate Std. Error t value Pr(>|t|)    
# k1 -5.558e-02  1.574e-02  -3.531 0.000574 ***
# k2 -3.348e-02  9.292e-03  -3.603 0.000448 ***
# k3 -2.700e-02  1.109e-02  -2.434 0.016318 *  
# A1  1.265e-02  7.481e-02   0.169 0.865979    
# A2 -1.035e-02  8.006e-02  -0.129 0.897290    
# A3  1.265e-05  7.859e-02   0.000 0.999872  
# NOTE: A not significant, use m2e

plot(resid(m2e)~ecklonia_decomp_21$time)
abline(0,0)
hist(resid(m2e))

m5e <- nls(rb ~ (exp(k[depth] * time)) + A[depth],
    start = list(k = c(-exp(-2.69769),-exp(-2.69769),-exp(-2.69769)),
                 A = c(0.09029,0.09029,0.09029)),
    data = ecklonia_decomp_21)
summary(m5e)
#m3e model with Residual fraction (A) does not fit the model. Using gnls to normalize for time|depth
m4e <- gnls(rb ~ exp(k * time),
           start = list(k = c(-exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
           params = list(k ~ depth),
           weights = varExp(form = ~time|depth),
           data = ecklonia_decomp_21)
plot(resid(m4e, type = "normalized") ~ time, data = ecklonia_decomp_21)
abline(0,0)
hist(resid(m4e, type = "normalized"))
summary(m4e) #ECKLONIA 10m (v 20m and 50m)
# Coefficients:
#   Value  Std.Error   t-value p-value
# k.(Intercept) -0.07733572 0.01321718 -5.851152  0.0000
# k.depth20      0.02808353 0.01398615  2.007953  0.0477
# k.depth50      0.02416140 0.01780968  1.356644  0.1783
intervals(m4e)

# change order of depth to get summary for 20m
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("20","50","10"))
m4e <- gnls(rb ~ exp(k * time),
            start = list(I = c(0.94357, 0.94357, 0.94357),
                          k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
summary(m4e) #ECKLONIA 20m 
confint(m4e)
# Coefficients:
#   Value   Std.Error    t-value p-value
# k.(Intercept) -0.04925215 0.004573675 -10.768615  0.0000
# k.depth50     -0.00392217 0.012783176  -0.306823  0.7597
# k.depth10     -0.02808156 0.013985740  -2.007871  0.0477
intervals(m4e)

# change order of depth to get summary for 50m
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("50","10","20"))
m4e <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
summary(m4e)
# Value  Std.Error   t-value p-value
# k.(Intercept) -0.05317432 0.01193696 -4.454594  0.0000
# k.depth10     -0.02416154 0.01780970 -1.356650  0.1783
# k.depth20      0.00392236 0.01278317  0.306838  0.7597
intervals(m4e)

funList <- list(~k1 - k2) #, ~k1 - k3, ~k2 - k3)
gnlht(m3s, funList)

###########################
m6e <- nls(rb ~ I * exp(k * time),
           start = list(k = c(-exp(-2.69769)),
                        I = c(0.94357)),
           data = ecklonia_decomp_21)
summary (m6e)

ecklonia_new <- data.frame(time = rep(seq(0, 50, by = 1)))
                          # depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new$fit <- predict(m6e, newdata =  ecklonia_new)
ggplot() + geom_point(data = ecklonia_decomp_21, aes(time, rb), na.rm=TRUE, position=position_dodge(width=2), size =2, alpha = 0.2) +
  geom_line(data = ecklonia_new, aes(time, fit), size=1, alpha=0.6) +
  #geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  #geom_errorbar(data = ecklonia_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  #geom_point(data = ecklonia_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  #scale_color_manual(values=c("#00BFC4", "lightskyblue", "royalblue4")) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title= "Ecklonia radiata") +
  theme_bw()+
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

confint(m6e)
###########################

# change order of depth to get summary for 50m
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("50","10","20"))
m4e <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
summary(m4e) #ECKLONIA 50m 
# Coefficients:
#   Value  Std.Error   t-value p-value
# k.(Intercept) -0.05317432 0.01193696 -4.454594  0.0000
# k.depth10     -0.02416154 0.01780970 -1.356650  0.1783
# k.depth20      0.00392236 0.01278317  0.306838  0.7597

# FIGURE 2.B ################################################

# reorder depth factor
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("10","20","50"))

# time must be continuous
decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$species <- factor(decomp_2021$species)
ecklonia_decomp_21 <- subset(decomp_2021, species == "Ecklonia radiata")
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth)

#calculate mean initial gC for E.radiata
mean(1.4028406, 1.4260579, 1.6029657) #[1] 1.402841
ecklonia_decomp_21$prC <- (ecklonia_decomp_21$gC/1.402841)
ecklonia_decomp_21$prC[is.na(ecklonia_decomp_21$prC)] <- 0
ecklonia_decomp_21$prC #column created for percent remaining carbon


# Start by calculating means and SE for each depth rb (ECKLONIA)
ecklonia_decomp_summary_21 <- ddply (ecklonia_decomp_21, c("depth", "time"), summarise,
                            N    = sum(!is.na(rb)),
                            mean = mean(rb, na.rm=TRUE),
                            sd   = sd(rb, na.rm = TRUE),
                            se   = sd / sqrt(N))
ecklonia_decomp_summary_21$depth <- factor(ecklonia_decomp_summary_21$depth, levels = c("10","20","50"))
ecklonia_decomp_summary_21
# make a fit line for decay estimated by model 4
ecklonia_new <- data.frame(time = rep(seq(0, 50, by = 1), 2),
                         depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new$depth <- factor(ecklonia_new$depth, levels = c("10", "20"))

ecklonia_new$fit <- predict(m4e, newdata =  ecklonia_new)
ecklonia_new

# Plot rb means +-se for each species with decomp_new$fit
ggplot() + geom_point(data = ecklonia_decomp_21, aes(time, rb, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=2), size =2, alpha = 0.2) +
  geom_line(data = ecklonia_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = ecklonia_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = ecklonia_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  scale_color_manual(values=c("#00BFC4", "lightskyblue", "royalblue4")) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title= "Ecklonia radiata") +
  theme_bw()+
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

# HALF OPEN PLOT STYLE
er_plot <- ggplot() + geom_point(data = ecklonia_decomp_21, aes(time, rb, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=2), size = 1.5, alpha = 0.4) +
  geom_line(data = ecklonia_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.0), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.25), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.75), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(1.0), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = ecklonia_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = ecklonia_decomp_summary_21, aes(time, mean, fill=depth, shape=depth),colour="black", na.rm = TRUE, pch = 21, size = 3) + # means added from summary stats
  scale_color_manual(values=c("#00BFC4", "dodgerblue", "darkgreen")) +
  scale_fill_manual(values =c("#00BFC4", "dodgerblue", "darkgreen")) +
  ylim(0, 1) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title="Ecklonia radiata") +
  theme_half_open() +
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
er_plot

# ANOVA
model <- aov(rb~time*depth, data=scyto_decomp_21)
summary(model)
TukeyHSD(model)

###################################################################
# MODEL1S - scytothalia %rb with depth ################################

decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$depth <- factor(decomp_2021$depth)
decomp_2021$species <- factor(decomp_2021$species)
decomp_2021$site <- factor(decomp_2021$site)
decomp_2021$gC <- (decomp_2021$C * decomp_2021$dw)/100
decomp_2021$D <- ((1-decomp_2021$rb)/decomp_2021$time)
scyto_decomp_21 <- subset(decomp_2021, species == "Scytothalia dorycarpa")
scyto_decomp <- read.csv ("scyto_decomp.csv")
scyto_decomp$depth <- factor (scyto_decomp$depth)
scyto_decomp$time <- scyto_decomp$days
scyto_decomp$rb <- scyto_decomp$biomass
scyto_decomp_21$C
scyto_decomp_21$gC
#calculate mean initial gC for S.dorycarpa
mean(1.6968091, 1.9180787, 1.7150520, 1.5681666, 1.6064144, 1.6093168, 1.8729970, 1.6686021, 1.7789071, 1.6720632) #[1] 1.696809
scyto_decomp_21$prC <- (scyto_decomp_21$gC/1.696809)
scyto_decomp_21$prC[is.na(scyto_decomp_21$prC)] <- 0
scyto_decomp_21$prC #column created for percent remaining carbon
# scyto_decomp <- read.csv ("2021 Marmion Kelp Litter Data - Scyto.csv")
# scyto_decomp$species <- as.factor(scyto_decomp$species)
# scyto_decomp$depth <- as.factor(scyto_decomp$depth)
#scyto_decomp_21$time <- numeric(scyto_decomp_21$time)
summary(scyto_decomp_21)
scyto_decomp_21

m1s <- nls(rb ~ SSasymp(time, A, I, logk), data = scyto_decomp_21)
summary(m1s)
# Formula: rb ~ SSasymp(time, A, I, logk)
# Parameters:
#       Estimate Std. Error t value Pr(>|t|)    
# A     0.10430    0.02319   4.498 1.94e-05 ***
# I     0.99883    0.03494  28.589  < 2e-16 ***
# logk -2.11538    0.09559 -22.129  < 2e-16 ***

# for gnls
#       Value  Std.Error   t-value p-value
# A     0.1043025 0.02318937   4.49786       0
# I     0.9988258 0.03493706  28.58929       0
# logk -2.1153834 0.09559140 -22.12943       0

# NEW MODEL with values from m1 for start for each depth
# write out the entire eqn: simplify by removing I since it is 1
m2s <- nls(rb ~ exp(k[depth] * time),
           start = list(k = c(-exp(-2.11538),-exp(-2.11538),-exp(-2.11538))),
           data = scyto_decomp_21)
summary(m2s)
# Formula: biomass ~ exp(k[depth] * time)
# Parameters:
#     Estimate Std. Error t value Pr(>|t|)    
# k1 -0.103400   0.007492  -13.80   <2e-16 ***
# k2 -0.069507   0.003857  -18.02   <2e-16 ***
# k3 -0.104805   0.007766  -13.49   <2e-16 ***
plot(resid(m2s)~scyto_decomp_21$time)
abline(0,0)
hist(resid(m2s))

m3s <- nls(rb ~ (exp(k[depth] * time)) + A[depth],
           start = list(k = c(-exp(-2.11479),-exp(-2.11479),-exp(-2.11479)),
                        A = c(0.10433, 0.10433 ,0.10433)),
           data = scyto_decomp)
summary(m3s)
# Estimate Std. Error t value Pr(>|t|)    
# k1 -9.885e-02  1.025e-02  -9.644   <2e-16 ***
# k2 -6.503e-02  6.114e-03 -10.636   <2e-16 ***
# k3 -9.501e-02  9.914e-03  -9.583   <2e-16 ***
# A1  1.690e-02  2.492e-02   0.678    0.499    
# A2  1.675e-02  2.340e-02   0.716    0.476    
# A3  1.269e-05  2.711e-02   0.000    1.000   
plot(resid(m3s)~scyto_decomp$days)
abline(0,0)
hist(resid(m3s))


m4s <- gnls(rb ~ exp(k * time) + A,
            start = list(k = c(-exp(-2.1153834),-exp(-2.1153834),-exp(-2.1153834)),
                         A = c(0.1043025,0.1043025,0.1043025)),
            params = list(k ~ depth,
                          A ~ depth),
            weights = varExp(form = ~ time|depth),
            data = scyto_decomp)
plot(resid(m4s, type = "normalized") ~ time, data = scyto_decomp)
abline(0,0)
hist(resid(m4s, type = "normalized"))
summary(m4s) #SCYTOTHALIA 10m (v 20m and 50m)
# Coefficients:
# Value   Std.Error    t-value p-value
# k.(Intercept) -0.09739559 0.008545654 -11.397090  0.0000
# k.depth20      0.03281864 0.010203039   3.216555  0.0016
# k.depth50      0.00238330 0.012028928   0.198131  0.8433
# A.(Intercept)  0.01090224 0.020732640   0.525849  0.5999
# A.depth20     -0.00362188 0.029443493  -0.123011  0.9023
# A.depth50     -0.01088955 0.030119449  -0.361546  0.7183
intervals(m4s)


# calculate without A but with gnls
m4s <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = scyto_decomp_21)
summary(m4s) #k.(Intercept) for 10m data
# Value   Std.Error    t-value p-value
# k.(Intercept) -0.10270762 0.007650771 -13.424480  0.0000
# k.depth20      0.03410146 0.008649161   3.942748  0.0002
# k.depth50     -0.00209691 0.009799462  -0.213982  0.8310
# change order of depth to get summary for 20m

scyto_decomp_21$depth <- factor(scyto_decomp_21$depth, levels = c("20","50","10"))
m4s <- gnls(rb ~ exp(k * time) + A,
            start = list(k = c(-exp(-2.1153834),-exp(-2.1153834),-exp(-2.1153834)),
                         A = c(0.1043025,0.1043025,0.1043025)),
            params = list(k ~ depth,
                          A ~ depth),
            weights = varExp(form = ~ time|depth),
            data = scyto_decomp)
summary(m4s) #scyto 20m 
# Coefficients:
# Value   Std.Error    t-value p-value
# k.(Intercept) -0.06860776 0.004034210 -17.006491   0e+00
# k.depth50     -0.03619677 0.007332803  -4.936280   0e+00
# k.depth10     -0.03410018 0.008649182  -3.942591   2e-04
intervals(m4s)

# change order of depth to get summary for 50m
scyto_decomp_21$depth <- factor(scyto_decomp_21$depth, levels = c("50","10","20"))
m4s <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = scyto_decomp_21)
summary(m4s) #scyto 50m 
# Coefficients:
# k.(Intercept) -0.10480453 0.006123329 -17.115614   0.000
# k.depth10      0.00209688 0.009799458   0.213980   0.831
# k.depth20      0.03619824 0.007332731   4.936529   0.000
intervals(m4s)

model <- aov(rb~time*depth, data=scyto_decomp_21)
summary(model)
TukeyHSD(model)



# FIGURE 2.B ################################################

# reorder depth factor
scyto_decomp_21$depth <- factor(scyto_decomp_21$depth, levels = c("10","20","50"))
scyto_decomp_21
# Start by calculating means and SE for each depth rb (scyto)
scyto_decomp_summary <- ddply (scyto_decomp_21, c("depth", "time"), summarise,
                                     N    = sum(!is.na(rb)),
                                     mean = mean(rb, na.rm=TRUE),
                                     sd   = sd(rb, na.rm = TRUE),
                                     se   = sd / sqrt(N))
#scyto_decomp_summary$time <- factor (scyto_decomp_summary$time)
scyto_decomp_summary
# make a fit line for decay estimated by model 4
scyto_new <- data.frame(time = rep(seq(0, 50, by = 1), 2),
                           depth = c(rep("10", 51), rep("20", 51)))
scyto_new$depth <- factor(scyto_new$depth)
#scyto_new$time <- factor(scyto_new$time)
scyto_new$fit <- predict(m2s, newdata =  scyto_new)
#decomp_new$time <- numeric(decomp_new$time)
scyto_new
# scyto_new %<>% mutate (low = fit * 0.75, high = fit * 1.25, index = seq(1,102))

#scyto_decomp_21$time <- as.factor(scyto_decomp_21$time)
# Plot rb means +-se for each species with decomp_new$fit
ggplot() + geom_point(data = scyto_decomp_21, aes(time, rb, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=3), size = 2, alpha = 0.2) +
  # geom_violin(data = scyto_decomp_21, aes(time, rb, fill=depth), position="dodge", alpha = 0.5) + 
  geom_line(data = scyto_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  # geom_ribbon(data = scyto_decomp, aes(time, fit, shading=depth, ymin = low, ymax = high), alpha = 0.1) +
  geom_errorbar(data = scyto_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = scyto_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  scale_color_manual(values=c("#F8766D", "violet", "darkorange")) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title= "Scytothalia dorycarpa") +
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
# save plot to pdf. saves in the current working directory. specify dimensions with width, height, and units.
ggsave("scyto_decomp_plot.pdf", width = 3, height = 3, units = "in")

sd_plot <- ggplot() + geom_point(data = scyto_decomp_21, aes(time, rb, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=2), size = 1.5, alpha = 0.4) +
  geom_line(data = scyto_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.0), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.25), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(0.75), linetype="dotted", alpha = 0.5) +
  geom_hline(yintercept=c(1.0), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = scyto_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = scyto_decomp_summary_21, aes(time, mean, fill=depth, shape=depth),colour="black", na.rm = TRUE, pch = 21, size = 3) + # means added from summary stats
  scale_color_manual(values=c("#00BFC4", "dodgerblue", "darkgreen")) +
  scale_fill_manual(values =c("#00BFC4", "dodgerblue", "darkgreen")) +
  ylim(0, 1) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title="Scytothalia dorycarpa") +
  theme_half_open() +
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
sd_plot

######################################################################
############## DECOMPOSITION CHEMISTRY ###############################
######################################################################

############################ GGPLOT FACET PLOTS ######################
decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$depth <- as.character(decomp_2021$depth)
decomp_2021$N <- as.numeric(decomp_2021$N)
decomp_2021$nitrogen <- decomp_2021$N
# decomp_2021$time <- as.character(decomp_2021$time)
decomp_2021$C.N <- as.numeric(decomp_2021$C.N)
decomp_2021$species <- factor(decomp_2021$species)
decomp_2021$time <- factor(decomp_2021$time)
scyto_decomp_21 <- subset(decomp_2021, species == "Scytothalia dorycarpa")
ecklonia_decomp_21 <- subset(decomp_2021, species == "Ecklonia radiata")
# summarise means +- SE for C, N, C/N, d15N


# % CARBON SUMMARY
# after summary is run for Ecklonia, switch the order of species to get summary for Scyto
lmC <- lm(C~time*species, data = decomp_2021)
summary(lmC)
mcheck(lmC)
decomp_2021$species <- c("Scytothalia dorycarpa", "Ecklonia radiata")
lmC <- lm(C~time*species, data = decomp_2021)
summary(lmC)
decomp_2021$species <- c("Ecklonia radiata", "Scytothalia dorycarpa")
# diff %C in 20 vs 10 and 50m
lmC <- lm(N~time*depth, data = ecklonia_decomp_21)
summary(lmC)
# no sd
chem_decomp_summary_21_C <- ddply (decomp_2021, c("species", "time"), summarise,
                                     N    = sum(!is.na(C)),
                                     mean = mean(C, na.rm=TRUE),
                                     sd   = sd(C, na.rm = TRUE),
                                     se   = sd / sqrt(N))

chem_decomp_summary_21_C
# species time  N     mean        sd        se
# 1      Ecklonia radiata    0  3 32.44333 0.7414400 0.4280706
# 2      Ecklonia radiata   14 49 30.59857 1.7777877 0.2539697
# 3      Ecklonia radiata   31 20 29.35600 1.5696879 0.3509929
# 4      Ecklonia radiata   50 10 30.10700 1.0021316 0.3169018
# 5 Scytothalia dorycarpa    0 10 31.56000 0.9696506 0.3066304
# 6 Scytothalia dorycarpa   14 54 30.87426 1.2068525 0.1642318
# 7 Scytothalia dorycarpa   31 22 31.05182 1.7420777 0.3714122
# 8 Scytothalia dorycarpa   50 11 33.15727 2.0524526 0.6188378

ecklonia_decomp_21_shelf <- ecklonia_decomp_21
ecklonia_decomp_21_shelf$location <- ecklonia_decomp_21_shelf$depth
levels(ecklonia_decomp_21_shelf$location) <- c("inshore", "inshore", "offshore")
ecklonia_decomp_21_shelf$location
ecklonia_decomp_21_shelf_summary <- ddply (ecklonia_decomp_21_shelf, c("location", "time"), summarise,
                                   N    = sum(!is.na(rb)),
                                   Q1   = quantile(rb, na.rm = TRUE,probs = 0.25),
                                   median=median(rb, na.rm = TRUE), 
                                   Q3   = quantile(rb, na.rm = TRUE, probs = 0.75),
                                   mean = mean(rb, na.rm=TRUE),
                                   sd   = sd(rb, na.rm = TRUE),
                                   se   = sd / sqrt(N))


ecklonia_decomp_21_shelf_summary
#### CALCS FOR % REMAINING carbon
# location time  N      mean         sd         se
# 1  inshore    0  3 1.0530688 0.07802549 0.04504803
# 2  inshore   14 35 0.3767720 0.27901686 0.04716246
# 3  inshore   31 21 0.2595491 0.22270473 0.04859816
# 4  inshore   50 15 0.0552703 0.10276635 0.02653416
# 5 offshore   14 18 0.4789360 0.36052170 0.08497578
##### CALCS BASED ON %C
# location time  N     mean       sd        se
# 1  inshore    0  3 32.44333 0.741440 0.4280706
# 2  inshore   14 32 30.78375 1.539689 0.2721812
# 3  inshore   31 20 29.35600 1.569688 0.3509929
# 4  inshore   50 10 30.10700 1.002132 0.3169018
# 5 offshore   14 17 30.25000 2.165823 0.5252891
### RB QUANTILES
# location time  N    Q1 median    Q3      mean         sd         se
# 1  inshore    0  3 0.950   0.98 1.035 0.9966667 0.08621678 0.04977728
# 2  inshore   14 35 0.135   0.44 0.550 0.3717143 0.26800445 0.04530102
# 3  inshore   31 21 0.050   0.23 0.480 0.2723810 0.23077922 0.05036015
# 4  inshore   50 15 0.000   0.03 0.055 0.0580000 0.11078938 0.02860569
# 5 offshore   14 18 0.195   0.53 0.655 0.4750000 0.34085101 0.08033935

scyto_decomp_21_shelf <- scyto_decomp_21
scyto_decomp_21_shelf$location <- factor(scyto_decomp_21_shelf$depth)
levels(scyto_decomp_21_shelf$location) <- c("inshore", "inshore", "offshore")
scyto_decomp_21_shelf$location
scyto_decomp_21_shelf_summary <- ddply (scyto_decomp_21_shelf, c("location", "time"), summarise,
                                           N    = sum(!is.na(C)),
                                           Q1   = quantile(C, na.rm = TRUE,probs = 0.25),
                                           median=median(C, na.rm = TRUE), 
                                           Q3   = quantile(C, na.rm = TRUE, probs = 0.75),
                                           mean = mean(C, na.rm=TRUE),
                                           sd   = sd(C, na.rm = TRUE),
                                           se   = sd / sqrt(N))


scyto_decomp_21_shelf_summary
# RB
# depth time  N     Q1 median     Q3      mean         sd         se
# 2  inshore   14 36 0.1900  0.270 0.3325 0.2858333 0.14059821 0.02343304
# 3  inshore   31 22 0.0900  0.130 0.1875 0.1477273 0.09879907 0.02106403
# 5 offshore   14 18 0.1675  0.225 0.2975 0.2305556 0.08495481 0.02002404


# Calculate interquartile ranges
ecklonia_decomp_21_shelf_14 <- subset(ecklonia_decomp_21_shelf, time == "14")
ecklonia_decomp_21_shelf_14in <- subset(ecklonia_decomp_21_shelf, location == "inshore")
ecklonia_decomp_21_shelf_14off <- subset(ecklonia_decomp_21_shelf, location == "offshore")
ecklonia_decomp_21_shelf_31 <- subset(ecklonia_decomp_21_shelf, time == "31")
ecklonia_decomp_21_shelf_50 <- subset(ecklonia_decomp_21_shelf, time == "50")
quantile(ecklonia_decomp_21_shelf_14in$prC, prob=c(.25,.5,.75), type=1)
quantile(ecklonia_decomp_21_shelf_14off$prC, prob=c(.25,.5,.75), type=1)
quantile(ecklonia_decomp_21_shelf_31$prC, prob=c(.25,.5,.75), type=1)
quantile(ecklonia_decomp_21_shelf_50$prC, prob=c(.25,.5,.75), type=1)

#                       25%        50%        75% 
# Ecklonia 14 (inshore) 0.03332409 0.20677090 0.51689015
# Ecklonia 31 (inshore) 0.04422454 0.20677090 0.45487678  
# Ecklonia 50 (inshore) 0.00000000 0.02642352 0.05354548 
# Ecklonia 14 (offshore)0.1778580 0.4803583 0.6257796 

ecklonia_decomp_21_shelf$time <- factor(ecklonia_decomp_21_shelf$time)
ggplot(ecklonia_decomp_21_shelf, aes(x=time, y=prC, fill=depth)) +
  geom_boxplot()
ggplot(ecklonia_decomp_21_shelf, aes(x=time, y=prC, fill=location)) +
  geom_boxplot()
# geom_violin(data = scyto_decomp_21, aes(time, rb, fill=depth), position="dodge", alpha = 0.5) + 
geom_line(data = scyto_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  # geom_ribbon(data = scyto_decomp, aes(time, fit, shading=depth, ymin = low, ymax = high), alpha = 0.1) +
  geom_errorbar(data = scyto_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = scyto_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  scale_color_manual(values=c("#F8766D", "violet", "darkorange")) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title= "Scytothalia dorycarpa") +
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
# save plot to pdf. saves in the current working directory. specify dimensions with width, height, and units.
ggsave("scyto_decomp_plot.pdf", width = 3, height = 3, units = "in")



insitu_export_summary <- ddply (ecklonia)
# % NITROGEN SUMMARY
# if error in evals, re read in data
lmN <- lm(nitrogen~species*time, data = decomp_2021)
summary(lmN)
mcheck(lmN)
decomp_2021$species <- factor(decomp_2021$species, levels = c("Scytothalia dorycarpa",
                                                              "Ecklonia radiata"))
lmN <- lm(nitrogen~species*time, data = decomp_2021)
summary(lmN)
decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata",
                                                              "Scytothalia dorycarpa"))

chem_decomp_summary_21_N <- ddply (decomp_2021, c("species", "time"), summarise,
                                      N    = sum(!is.na(nitrogen)),
                                      mean = mean(nitrogen, na.rm=TRUE),
                                      sd   = sd(nitrogen, na.rm = TRUE),
                                      se   = sd / sqrt(N))

chem_decomp_summary_21_N
# species time  N      mean         sd         se
# 1      Ecklonia radiata    0  3 0.4400000 0.17578396 0.10148892
# 2      Ecklonia radiata   14 49 0.7689796 0.37815476 0.05402211
# 3      Ecklonia radiata   31 20 0.7175000 0.26643603 0.05957691
# 4      Ecklonia radiata   50 10 1.0040000 0.20457001 0.06469072
# 5 Scytothalia dorycarpa    0  8 0.1600000 0.09913915 0.03505098
# 6 Scytothalia dorycarpa   14 53 0.2413208 0.12863242 0.01766902
# 7 Scytothalia dorycarpa   31 22 0.1681818 0.07979953 0.01701332
# 8 Scytothalia dorycarpa   50 11 0.1836364 0.04522670 0.01363636

# CARBON/NITROGEN SUMMARY
# lm summary for C/N in Ecklonia (timepoint 0 vs 14, 31, 50)
lmCN <- lm(C.N~species*time, data = decomp_2021)
summary(lmCN)
mcheck(lmCN)

# lm summary for d15N in Scyto (timepoint 0 vs 14, 31, 50) - reordered species
decomp_2021$species <- factor(decomp_2021$species, levels = c("Scytothalia dorycarpa",
                                                              "Ecklonia radiata"))
lmd15N <- lm(d15N~species*time, data = decomp_2021)
summary(lmd15N)
# reorder time to look at T2 vs T3
decomp_2021$time <- factor(decomp_2021$time, levels = c("14", "0", "31", "50"))
summary(lmd15N)
# reset order to original configuration (alphabetical and ascending)
decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata",
                                                              "Scytothalia dorycarpa"))
decomp_2021$time <- factor(decomp_2021$time, levels = c("0", "14", "31", "50"))
chem_decomp_summary_21_C.N <- ddply (decomp_2021, c("species", "time"), summarise,
                                   N    = sum(!is.na(C.N)),
                                   mean = mean(C.N, na.rm=TRUE),
                                   sd   = sd(C.N, na.rm = TRUE),
                                   se   = sd / sqrt(N))

chem_decomp_summary_21_C.N
# species time  N      mean         sd        se
# 1      Ecklonia radiata    0  3  85.59333  45.544256 26.294988
# 2      Ecklonia radiata   14 49  54.00612  39.004150  5.572021
# 3      Ecklonia radiata   31 20  49.55000  28.817291  6.443742
# 4      Ecklonia radiata   50 10  31.51100   8.505168  2.689570
# 5 Scytothalia dorycarpa    0  8 324.32875 252.697155 89.341936
# 6 Scytothalia dorycarpa   14 53 196.84264 186.676074 25.641931
# 7 Scytothalia dorycarpa   31 22 240.02318 138.699733 29.570882
# 8 Scytothalia dorycarpa   50 11 194.13909  61.861856 18.652052

# NITROGEN ISOTOPE d15N SUMMARY
# lm summary for d15N in Ecklonia (timepoint 0 vs 14, 31, 50)
lmd15N <- lm(d15N~species*time, data = decomp_2021)
summary(lmd15N)
mcheck(lmd15N)

# lm summary for d15N in Scyto (timepoint 0 vs 14, 31, 50) - reordered species
decomp_2021$species <- factor(decomp_2021$species, levels = c("Scytothalia dorycarpa",
                                                              "Ecklonia radiata"))
lmd15N <- lm(d15N~species*time, data = decomp_2021)
summary(lmd15N)
# reorder time to look at T2 vs T3
decomp_2021$time <- factor(decomp_2021$time, levels = c("14", "0", "31", "50"))
summary(lmd15N)
# reset order to original configuration (alphabetical and ascending)
decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata",
                                                              "Scytothalia dorycarpa"))
decomp_2021$time <- factor(decomp_2021$time, levels = c("0", "14", "31", "50"))

# create a summary of mean values for fig and table
chem_decomp_summary_21_d15N <- ddply (decomp_2021, c("species", "time"), summarise,
                                     N    = sum(!is.na(d15N)),
                                     mean = mean(d15N, na.rm=TRUE),
                                     sd   = sd(d15N, na.rm = TRUE),
                                     se   = sd / sqrt(N))

chem_decomp_summary_21_d15N
# species time  N     mean        sd         se
# 1      Ecklonia radiata    0  3 3.510000 0.4703190 0.27153882
# 2      Ecklonia radiata   14 49 4.668163 0.5018701 0.07169573
# 3      Ecklonia radiata   31 20 4.625000 0.3781186 0.08454989
# 4      Ecklonia radiata   50 10 4.602000 0.3807828 0.12041410
# 5 Scytothalia dorycarpa    0 10 2.685000 0.6714536 0.21233229
# 6 Scytothalia dorycarpa   14 54 3.405185 0.3387309 0.04609544
# 7 Scytothalia dorycarpa   31 22 3.000000 0.4585069 0.09775400
# 8 Scytothalia dorycarpa   50 11 3.136364 0.4694097 0.14153235

chem_decomp_summary_21_d13C <- ddply (decomp_2021, c("species", "time"), summarise,
                                      N    = sum(!is.na(d13C)),
                                      mean = mean(d13C, na.rm=TRUE),
                                      sd   = sd(d13C, na.rm = TRUE),
                                      se   = sd / sqrt(N))

chem_decomp_summary_21_d13C

# DECAY RATE
lmD <- lm(D~species*time, data = decay_decomp_21)
summary(lmD)
mcheck(lmN)
chem_decomp_summary_21_D <- ddply (decay_decomp_21, c("species", "time"), summarise,
                                   N    = sum(!is.na(D)),
                                   mean = mean(D, na.rm=TRUE),
                                   sd   = sd(D, na.rm = TRUE),
                                   se   = sd / sqrt(N))



chem_decomp_summary_21_D
# 1      Ecklonia radiata   14 53 0.04308625 0.021114815 0.0029003430
# 2      Ecklonia radiata   31 21 0.02379416 0.007444491 0.0016245211
# 3      Ecklonia radiata   50 15 0.01904000 0.002215788 0.0005721139
# 4 Scytothalia dorycarpa   14 54 0.05304233 0.009052307 0.0012318630
# 5 Scytothalia dorycarpa   31 22 0.02781525 0.003187067 0.0006794849
# 6 Scytothalia dorycarpa   50 12 0.01870000 0.001283461 0.0003705033

decomp_2021$time <- factor(decomp_2021$time)
decomp_2021 %>% 
  ggplot(aes(time, C, fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  ylim(27,37.5) +
  theme_bw() +
  labs(x="Detrital age (days)", y= "% Carbon") +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

carbon2021 <- chem_decomp_summary_21_C %>% 
  ggplot (aes(time, mean, group = species, color = species), colour="black", na.rm = TRUE, pch = 21) +
  geom_line () +
  geom_errorbar(data = chem_decomp_summary_21_C, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  geom_point(aes(time, mean, group = species, fill = species), colour="black", na.rm = TRUE, pch = 21, size = 3, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#F8766D","#00BFC4")) +
  labs(x="Detrital age (days)", y= "% Carbon (of DW)") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
carbon2021

mcheck(mC)
  
decomp_2021 %>% 
  ggplot(aes(time, D, fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  #ylim(-0.02,0.15) +
  theme_bw() +
  labs(x="Detrital age (days)", y= "Decomposition (%/day)") +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

chem_decomp_summary_21_D %>% 
  ggplot(aes(time, mean, group = species, color = species, shape = species), na.rm = TRUE) +
  geom_line () +
  geom_point(size = 3, position=position_dodge(width=.1)) +
  geom_errorbar(data = chem_decomp_summary_21_D, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  labs(x="Detrital age (days)", y= "Decomposition (%biomass/day)") +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

chem_decomp_summary_21_D %>% 
  ggplot(aes(time, mean, group = species, color = species, shape = species), na.rm = TRUE) +
  geom_line () +
  geom_point(size = 3, position=position_dodge(width=.1)) +
  geom_errorbar(data = chem_decomp_summary_21_D, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  labs(x="Detrital age (days)", y= "Decomposition (%biomass/day)") +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

decomp_2021 %>% 
  ggplot(aes(time, N, fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  ylim(0, 1.5) +
  theme_bw() +
  labs(x="Detrital age (days)", y= "% Nitrogen") +
  theme(legend.position = 'bottom', legend.direction = "horizontal")
  
chem_decomp_summary_21_N %>% 
  ggplot(aes(time, mean, group = species, color = species, shape = species), na.rm = TRUE) +
  geom_line () +
  geom_point(size = 3, position=position_dodge(width=.1)) +
  geom_errorbar(data = chem_decomp_summary_21_N, aes(time, mean, ymin = mean - sd, ymax = mean + sd), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  labs(x="Detrital age (days)", y= "% Nitrogen") +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = "horizontal") 
nitrogen2021 <- chem_decomp_summary_21_N %>% 
  ggplot (aes(time, mean, group = species, color = species), colour="black", na.rm = TRUE, pch = 21) +
  geom_line () +
  geom_errorbar(data = chem_decomp_summary_21_N, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  geom_point(aes(time, mean, group = species, fill = species), colour="black", na.rm = TRUE, pch = 21, size = 3, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#F8766D","#00BFC4")) +
  labs(x="Detrital age (days)", y= "% Nitrogen (of DW)") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(legend.text = element_text(face = "italic"))
nitrogen2021
ggsave("Nplot", nitrogen2021, width = 3, height = 3)
#combine C and N plots into one plot
combined_plots <- (carbon2021 + nitrogen2021) / theme(legend.text = element_text(face = "italic"))

combined_plots

decomp_2021 %>% 
  ggplot(aes(time, C.N, fill = species)) +
  geom_boxplot() +
  ylim(0,750) +
  facet_wrap(~species) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  ylim(0,200)+
  theme_bw() +
  labs(x="Detrital age (days)", y= "Carbon-nitrogen ratio") +
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme_half_open()

chem_decomp_summary_21_C.N %>% 
  ggplot(aes(time, mean, group = species, color = species), na.rm = TRUE) +
  geom_line () +
  geom_point(size = 3, position=position_dodge(width=.1)) +
  geom_errorbar(data = chem_decomp_summary_21_C.N, aes(time, mean, ymin = mean - sd, ymax = mean + sd), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#F8766D","#00BFC4")) +
  labs(x="Detrital age (days)", y= "Carbon-nitrogen ratio") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal")
chem_decomp_summary_21_C.N %>% 
  ggplot (aes(time, mean, group = species, color = species), colour="black", na.rm = TRUE, pch = 21) +
  geom_line () +
  geom_errorbar(data = chem_decomp_summary_21_C.N, aes(time, mean, ymin = mean - sd, ymax = mean + sd), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  geom_point(aes(time, mean, group = species, fill = species), colour="black", na.rm = TRUE, pch = 21, size = 3, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#F8766D","#00BFC4")) +
  labs(x="Detrital age (days)", y= "Carbon-nitrogen ratio") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))

decomp_2021 %>% 
  ggplot(aes(time, d15N, fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  ylim(2,6)+
  theme_bw() +
  labs(x="Detrital age (days)", y= "δ15N (µm/g)") +
  theme(legend.position = 'bottom', legend.direction = "horizontal")
# drop in percent N due to leaching of organic vol comp? (lit search)

chem_decomp_summary_21_d15N %>% 
  ggplot(aes(time, mean, group = species, color = species, shape = species), na.rm = TRUE) +
  geom_line () +
  geom_errorbar(data = chem_decomp_summary_21_d15N, aes(time, mean, ymin = mean - sd, ymax = mean + sd), na.rm = TRUE, width = 0.1, position=position_dodge(width=.1)) +
  geom_point(aes(time, mean, group = species, color = species, fill = species), colour="black", pch = 21, size = 3, position=position_dodge(width=.1)) +
  scale_color_manual(values=c("#F8766D","#00BFC4")) +
  labs(x="Detrital age (days)", y= "δ15N (µm/g)") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))

lmdd15N <- lm(d15N~species*time, data = decomp_2021)
mcheck(lmdd15N)
summary(lmdd15N)

decomp_2021 %>% 
  ggplot(aes(time, d13C,  fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  theme_bw()

chem_decomp_summary_21_d13C %>% 
  ggplot(aes(time, mean, group = species, color = species, shape = species), na.rm = TRUE) +
  geom_line () +
  geom_point(size = 3) +
  geom_errorbar(data = chem_decomp_summary_21_d13C, aes(time, mean, ymin = mean - sd, ymax = mean + sd), na.rm = TRUE, width = 0.1) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  labs(x="Detrital age (days)", y= "δ13C (µm/g)") +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = "horizontal")

lmd13C <- lm(d13C~species*time, data = decomp_2021)
mcheck(lmd13C)
summary(lmd13C)
decomp_2021$species <- factor(decomp_2021$species, levels = c("Scytothalia dorycarpa",
                                                              "Ecklonia radiata"))
lmd13C <- lm(d13C~species*time, data = decomp_2021)
mcheck(lmd13C)
summary(lmd13C)
decomp_2021$species <- factor(decomp_2021$species, levels = c("Ecklonia radiata",
                                                              "Scytothalia dorycarpa"))



##########

library(lme4)

# DW by DEPTH (ECKLONIA)
df_Er <- decomp_2021 %>% 
  filter(species =="Ecklonia radiata") %>% 
  group_by(depth) 


l1 <- lm(C.N ~ time, data = df_Er)

anova(l1)
TukeyHSD(l1)

#            Df Sum Sq Mean Sq F value Pr(>F)
#depth        1  0.988  0.9881   0.913  0.352
#Residuals   18 19.492  1.0829 

# visualise residuals
plot(resid(m3)~decomp_2021$time)
abline(0,0)
hist(resid(m3))

# Model as data frame
ecklonia_new <- data.frame(time = rep(seq(0, 50, by = 1), 2),
                           depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new$depth <- factor(ecklonia_new$depth)
ecklonia_new$fit <- predict(m3e, newdata =  ecklonia_new)

ecklonia_new$time <- numeric(ecklonia_new$time)
ggplot() + geom_point(data = ecklonia_decomp_21, aes(time, rb, color=depth)) + 
  geom_line(data = ecklonia_new, aes(time, fit, color=depth)) +
  theme_bw()

# STOP
# STOP
# STOP
# STOP
# STOP
# STOP
# STOP
# STOP
# STOP

#### ANOVA WITH MARTY ###########################################
aov1=aov(rb~factor(time)*species, data=decomp_2021)
summary(aov1)
TukeyHSD(aov1)



mcheck <- function (obj, ... ) {
  rs<-resid(obj)
  fv<-fitted(obj)
  par(mfrow=c(1,2))
  require(car)
  plot(fv,rs,xlab="Fitted values",ylab="Residuals")
  abline(h=0, lty=2)
  lines(smooth.spline(fv, rs), col = "red")
  qqPlot(rs,xlab="Normal scores",ylab="Ordered residuals")
  par(mfrow=c(1,1))
  invisible(NULL)}
mcheck(m4e)


table(ecklonia_decomp$depth, ecklonia_decomp$site)


aov1=aov(dw~factor(time)*site, data=ecklonia_decomp)
summary(aov1)
TukeyHSD(aov1)

table(ecklonia_decomp$depth, ecklonia_decomp$site)













m2 <- nls(rb ~ exp(k[species] * time) + A[species],
          start = list(k = c(-exp(-2.33516),-exp(-2.33516)),
                       A = c(0.11, 0.11)),
          data = decomp_2021)
summary(m3)
m1e <- nls(biomass ~ exp(k[depth] * time) + A[depth],
           start = list(k = c(-exp(-3.55552),-exp(-3.55552),-exp(-3.55552)),
                        A = c(-0.16768,-0.16768,-0.16768)),
           data = ecklonia_decomp)
summary(m3e)


########## (%RB) - ECKLONIA SUMMARY DATA ########################
#################################################################

ecklonia_decomp_21 <- subset (decomp_2021, species == "Ecklonia radiata")
ecklonia_decomp_21$time <- factor(ecklonia_decomp_21$time)
# Run the functions mean, sd, and se on the value of "percent.biomass.remaining" for each group, 
# broken down by depth + time

ggplot() + geom_point(data = decomp_202)

ecklonia_decomp_summary_21 <- ddply (ecklonia_decomp_21, c("depth", "time"), summarise,
                                     N    = sum(!is.na(rb)),
                                     mean = mean(rb, na.rm=TRUE),
                                     sd   = sd(rb, na.rm = TRUE),
                                     se   = sd / sqrt(N))
ecklonia_decomp_summary_21
# depth time  N       mean         sd         se
# 1    10    0  3 1.00333333 0.03785939 0.02185813
# 2    10   14 19 0.32157895 0.23596759 0.05413468
# 3    10   31  4 0.19500000 0.23014488 0.11507244
# 4    10   50  3 0.15000000 0.25980762 0.15000000
# 5    20   14 16 0.46125000 0.31129568 0.07782392
# 6    20   31 17 0.55000000 0.67345193 0.16333608
# 7    20   50 12 0.03666667 0.03498918 0.01010051
# 8    50   14 18 0.49333333 0.35418506 0.08348222

ecklonia_decomp_summary_21$time
#[1] 0  14 31 50 14 31 50 14

ecklonia_decomp_summary_21$mean
# [1] 0.99666667 0.30894737 0.18750000 0.14666667 0.44625000 0.29235294 0.03583333 0.47500000

# summarising data using summarySE fn w 95% ci
# need to rework: summarySE(ecklonia_decomp_21, measurevar="biomass", groupvars=c("depth","time"))
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

# DATA PREP
ecklonia_decomp_summary_21$depth <- as.factor(ecklonia_decomp_summary_21$depth)
# ecklonia_decomp_summary_21$time <- as.integer(ecklonia_decomp_summary_21$time)



########## (%RB) - SCYTOTHALIA SUMMARY DATA #####################
#################################################################

scytothalia_decomp_21 <- subset (decomp_2021, species == "Scytothalia dorycarpa")
scytothalia_decomp_21

# Run the functions mean, sd, and se on the value of "percent.biomass.remaining" for each group, 
# broken down by depth + time
# need to address NAs 

scytothalia_decomp_summary_21 <- ddply (scytothalia_decomp_21, c("depth", "time"), summarise,
                                     N    = sum(!is.na(rb)),
                                     mean = mean(rb, na.rm=TRUE),
                                     sd   = sd(rb, na.rm = TRUE),
                                     se   = sd / sqrt(N))
scytothalia_decomp_summary_21

# depth time  N      mean         sd         se
# 1    10    0 10 1.0010000 0.12224293 0.03865661
# 2    10   14 18 0.2222222 0.07689225 0.01812368
# 3    10   31  4 0.0950000 0.05802298 0.02901149
# 4    10   50  2 0.0900000 0.12727922 0.09000000
# 5    20   14 18 0.3316667 0.16165504 0.03810246
# 6    20   31 18 0.1538889 0.10024317 0.02362754
# 7    20   50 10 0.0690000 0.05384133 0.01702612
# 8    50   14 18 0.2250000 0.08197919 0.01932268

scytothalia_decomp_summary_21$time
#[1] 0  14 31 50 14 31 50 14

scytothalia_decomp_summary_21$mean
# [1] 1.0010000 0.2222222 0.0950000 0.0900000 0.3316667 0.1538889 0.0690000 0.2250000

# summarising data using summarySE fn w 95% ci
# need to rework: summarySE(scytothalia_decomp_21, measurevar="biomass", groupvars=c("depth","time"))
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

# DATA PREP
scytothalia_decomp_summary_21$depth <- factor(scytothalia_decomp_summary_21$depth)
# scytothalia_decomp_summary_21$time <- as.integer(scytothalia_decomp_summary_21$time)


####### CHECK FOR INTERACTIONS BETWEEN SITE AND DEPTH ###########
#################################################################
# Visualise the data

decomp_2021$time <- factor(decomp_2021$time)
decomp_2021 %>% 
  select(species, site, depth, time, rb) %>% 
  # filter(species == "Ecklonia radiata") %>% 
  na.omit() %>% 
  ggplot(aes(time, rb,
             color = depth))+
  geom_boxplot(size = 1, alpha = 0.5) + # alpha is transparency +
  # geom_hline(yintercept=c(0,10), linetype="dotted", alpha = 0.5) +
  facet_wrap(~species) + 
  # ylim(-2.6,11) +
  labs(title = "MARMION IN SITU KELP DECOMPOSITION EXPERIMENT - DEPTH EFFECTS?")+
  theme_bw()

# (?) Scytothalia samples collected at DEEP 3 may have substantially less remaining biomass (DW) than DEEP1 or DEEP2 
# (?) Syctothalia samples collected at MID sites - possibly a slower decay than 10 & 50m

# 2 WAY ANOVA: ECKLONIA REMAINING BIOMASS OVER TIME (FACTORS= TIME (x4) AND DEPTH(x3) + NESTED SITE(x3 ant each depth))


# Visualise homoscedasticity
par(mfrow=c(2,2))
plot(aov1e)
par(mfrow=c(1,1))

a1 <- lm(-log(rb)~factor(time)*species, data = decomp_2021)
summary(a1)

mcheck <- function (obj, ... ) {
  rs<-resid(obj)
  fv<-fitted(obj)
  par(mfrow=c(1,2))
  require(car)
  plot(fv,rs,xlab="Fitted values",ylab="Residuals")
  abline(h=0, lty=2)
  lines(smooth.spline(fv, rs), col = "red")
  qqPlot(rs,xlab="Normal scores",ylab="Ordered residuals")
  par(mfrow=c(1,1))
  invisible(NULL)}
mcheck(m4)

Anova(lm(rb~time*depth, data = ecklonia_decomp_21, contrasts=list(time="contr.sum")), type="2")
Anova(aov1e, type = "II")
# Anova Table (Type II tests)
# Response: rb
#                     Sum Sq Df F value    Pr(>F)    
# factor(time)       2.5591  3 13.1713 3.967e-07 ***
# depth              0.2387  2  1.8430    0.1647    
# factor(time):depth 0.1169  2  0.9026    0.4094    
# Residuals          5.4402 84 


# Post-hoc test
TukeyHSD(aov1)
# $`factor(time)`
# diff        lwr         upr     p adj
# 14-0  -0.5898742 -0.9857565 -0.19399190 0.0010686
# 31-0  -0.7242857 -1.1360095 -0.31256190 0.0000827
# 50-0  -0.9386667 -1.3605580 -0.51677530 0.0000006
# 31-14 -0.1344115 -0.3064156  0.03759265 0.1788407
# 50-14 -0.3487925 -0.5438854 -0.15369952 0.0000620
# 50-31 -0.2143810 -0.4398914  0.01112947 0.0685025
# NOTES: Ecklonia rb is not significantly different between 14 and 31 or 31 and 50days

aov2e <- lm(rb~factor(time)*site, data = ecklonia_decomp_21, contrasts=list(time="contr.sum"), type="2")
summary(aov2e)

Anova(aov2e, type = "II")
# Anova Table (Type II tests)
# Response: rb
#                    Sum Sq Df F value   Pr(>F)   
# factor(time)      0.9373  2  7.1534 0.001434 **
# site              0.7317  8  1.3961 0.212167   
# factor(time):site 0.1508  5  0.4603 0.804520   
# Residuals         4.9134 75 

TukeyHSD(aov2e)
# NOTE: NO significant difference within rb of sites at the same depth. Supported by initial visualisation boxplot

# 2 WAY ANOVA: SCYTOTHALIA PERCENT REMAINING BIOMASS OVER TIME (FACTORS= TIME (x4) AND DEPTH(x3) + NESTED SITE(x3 ant each depth))
aov1s <- aov(rb~factor(time)*depth, data = scytothalia_decomp_21)
Anova(aov1s, type="II")
# Anova Table (Type II tests)
# Response: rb
# Sum Sq Df  F value    Pr(>F)    
# factor(time)       5.1657  3 159.2886 < 2.2e-16 ***
# depth              0.1359  2   6.2843  0.002788 ** 
# factor(time):depth 0.0286  2   1.3212  0.271930    
# Residuals          0.9729 90   

TukeyHSD(aov1s)
# $`factor(time)`
#           diff        lwr         upr     p adj
# 14-0  -0.74137037 -0.8371537 -0.64558702 0.0000000
# 31-0  -0.85781818 -0.9639292 -0.75170714 0.0000000
# 50-0  -0.92850000 -1.0476290 -0.80937096 0.0000000
# 31-14 -0.11644781 -0.1868191 -0.04607652 0.0002218
# 50-14 -0.18712963 -0.2759232 -0.09833609 0.0000020
# 50-31 -0.07068182 -0.1705287  0.02916504 0.2557502
# NOTES: no diff between 31 and 50 days for %rm~time for Scytothalia

# $depth
#           diff          lwr         upr     p adj
# 20-10  0.05606022 -0.001225757 0.113346194 0.0564743
# 50-10 -0.01018667 -0.084019185 0.063645844 0.9421891
# 50-20 -0.06624689 -0.136666969 0.004173191 0.0696703
# NOTES: no actual effect

aov2s <- aov(rb~factor(time)*site, data = scytothalia_decomp_21)
Anova(aov2s, type = "II")
# Anova Table (Type II tests)
#                   Df Sum Sq Mean Sq F value Pr(>F)    
# factor(time)       3  6.068  2.0228 192.801 <2e-16 ***
# site               8  0.217  0.0272   2.588 0.0143 *  
# factor(time):site  5  0.070  0.0141   1.340 0.2558    
# Residuals         81  0.850  0.0105 



TukeyHSD(aov2s)
# MID 2-DEEP 3         0.159320658  0.001074769 0.31756655 0.0470743
# NOTE: significant effect of site on rb above is attributed to differences between mid2 and deep3

# 14:DEEP 3-14:DEEP 1       -0.105000000 -0.343918449  0.133918449 0.9989216
# 14:DEEP 3-14:DEEP 2       -0.111666667 -0.350585115  0.127251782 0.9967710
# NOTE: no variation in rb between deep3 and deep2 or deep1 like boxplot suggested


aov1 <-aov(rb~factor(time)*species, data = decomp_2021)
summary(aov1)
TukeyHSD(aov1)




#################################################################
#   LINE GRAPH   ################################################
#################################################################

library(ggplot2)

# DATA PREP
ecklonia_decomp_summary_21$time <- as.factor(ecklonia_decomp_summary_21$time)
ecklonia_decomp_summary_21$depth <- as.factor(ecklonia_decomp_summary_21$depth)

pd <- position_dodge(0.1) # move them .05 to the left and right 

ggplot(ecklonia_decomp_summary_21, aes(x=time, y=mean, colour=depth, group=depth)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
  #geom_line(position=pd) + # lines between means
  #geom_line(data = ecklonia_new, aes(time, fit, color=depth))+ # m2e model to zero biomass
  geom_line(data = ecklonia_new, aes(time, fit, color=depth)) +
  scale_color_manual (values=c ('#ADD8E6', '#528B8B', '#00688B')) +
  geom_point(position=pd, size=2, shape=21, fill="grey") +
  geom_smooth(method = "gam", formula = y ~ poly(x, 2)) +
  xlab("Time(days)") +
  ylab("Percent remaining biomass (E.radiata)") +
  ggtitle("E. radiata in situ decay rate (wet weight)") +
  scale_colour_hue(name="Depth",    # Legend label, use darker colors
                   breaks=c("10", "20", "50"),
                   labels=c("10m", "20m", "50m"),
                   l=60) +                    # Use darker colors, lightness=40
  expand_limits(y=0) +                        # Expand y range
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(.95,.7))            # Position legend in bottom right

# geom_line(data=ecklonia_decomp_21_tmp, aes(x=time, y=pred))

#################################################################
########## PREDICTIVE MODEL FOR EXP DECAY #######################
#################################################################

# SSasymp is a self starting function
# C is the area above the asymptote
# A is the area below the asymptote
# I should be close to 1


#################################################################
#   chem analysis  ########################
#################################################################
decomp_2021$N <- as.numeric(decomp_2021$N)
decomp_2021$C.N <- as.numeric(decomp_2021$C.N)
decomp_2021$time <- as.factor(decomp_2021$time)

# PERCENT CARBON
decomp_2021 %>% 
  ggplot(aes(time, C, fill=species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  ylim(27,37) +
  theme_bw()
aovC <- aov(C~factor(time)*species, data = decomp_2021)
Anova(aovC, type="II")
TukeyHSD(aovC)
#Scytothalia dorycarpa-Ecklonia radiata 0.9410749 0.4880823 1.394068 6.36e-05
aovCE <- aov(rb~factor(time)*depth, data= scytothalia_decomp_21)
Anova(aovCE, type="II")
summary(aovCE)
# Response: C
# Sum Sq Df F value  Pr(>F)  
# factor(time)        31.547  3  3.9140 0.01188 *
# depth                3.175  2  0.5909 0.55643  
# factor(time):depth   6.670  2  1.2414 0.29494  
# Residuals          198.813 74 
TukeyHSD(aovCE)
aovCS <- aov(C~factor(time)*depth, data= scytothalia_decomp_21)
Anova(aovCS, type="II")
# Response: C
# Sum Sq Df F value  Pr(>F)  
# factor(time)        31.547  3  3.9140 0.01188 *
# depth                3.175  2  0.5909 0.55643  
# factor(time):depth   6.670  2  1.2414 0.29494  
# Residuals          198.813 74 
TukeyHSD(aovCE)
# Response: C
# Sum Sq Df F value    Pr(>F)    
# factor(time)        40.339  3  6.4172 0.0005513 ***
# depth                2.242  2  0.5351 0.5875089    
# factor(time):depth   2.785  2  0.6646 0.5170231    
# Residuals          186.486 89 

# PERCENT NITROGEN
decomp_2021 %>% 
  ggplot(aes(time, N, fill=species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  ylim(0, 1.5) +
  theme_bw()
aovN <- aov(N~factor(time)*species, data = decomp_2021)
Anova(aovN, type="II")
# Response: N
# Sum Sq  Df  F value  Pr(>F)    
# factor(time)          0.4936   3   2.8393 0.03960 *  
# species              13.3932   1 231.1383 < 2e-16 ***
# factor(time):species  0.5532   3   3.1824 0.02541 *  
# Residuals             9.7347 168   
TukeyHSD(aovN)
#Scytothalia dorycarpa-Ecklonia radiata -0.5502265 -0.6220356 -0.4784175     0
aovNE <- aov(N~factor(time)*depth, data= ecklonia_decomp_21)
Anova(aovNE, type="II")
TukeyHSD(aovNE)

aovNS <- aov(N~factor(time)*depth, data= scytothalia_decomp_21)
Anova(aovNS, type="II")
TukeyHSD(aovNS)

# CARBON/NITROGEN
decomp_2021 %>% 
  ggplot(aes(time, C.N, fill=species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  ylim(0,750) +
  theme_bw()
aovC.N <- aov(C.N~factor(time)*species, data = decomp_2021)
Anova(aovC.N, type="II")
# Response: C.N
# Sum Sq  Df F value    Pr(>F)    
# factor(time)          108259   3  2.1691   0.09356 .  
# species              1131738   1 68.0271 4.465e-14 ***
#   factor(time):species   30685   3  0.6148   0.60631    
# Residuals            2794945 168  
aovCNE <- aov(C.N~factor(time)*depth, data= ecklonia_decomp_21)
Anova(aovCNE, type="II")
TukeyHSD(aovC.NE)

#################################################################
#   smooth fit curve with SE envelopes   ########################
#################################################################

# plot(ecklonia_decomp_21$time,ecklonia_decomp_21$biomass, colour=ecklonia_decomp_21$depth, group=ecklonia_decomp_21$depth)

# https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
# make a new data frame for 20m 
ecklonia_mid <- subset(ecklonia_decomp, depth=="20")
ecklonia_mid$depth <- as.factor(ecklonia_mid$depth)
ecklonia_mid

ggplot(ecklonia_mid, aes(time, biomass)) + geom_point() +
  geom_smooth(method = "gam", formula = y ~ poly(x, 2)) 


#################################################################
#   2 WAY ANOVA   ###############################################
#################################################################

# DATA PREP
# Show a random sample
set.seed(1234)
dplyr::sample_n(ecklonia_decomp_21, 10)

# Check the structure
str(ecklonia_decomp_21)
# 'data.frame':	143 obs. of  4 variables:
# $ site   : chr  "shallow1" "shallow1" "shallow1" "shallow1" ...
# $ depth  : int  10 10 10 10 10 10 20 20 20 20 ...
# $ time   : int  0 0 0 0 0 0 0 0 0 0 ...
# $ biomass: num  1 1 1 1 1 1 1 1 1 1 ...

# Convert depth as a factor and recode the levels as 10m, 20m, 50m
ecklonia_decomp_21$depth <- as.factor(ecklonia_decomp_21$depth)
#below isn't working.... labels?
ecklonia_decomp_21$depth <- as.factor(ecklonia_decomp_21$depth,
                                      levels=c(10, 20, 50),
                                      labels=c("10m", "20m", "50m"))
ecklonia_decomp_21$time <- as.integer(ecklonia_decomp_21$time)
head(ecklonia_decomp_21)

# generate frequency tables
table(ecklonia_decomp_21$time, ecklonia_decomp_21$depth)
#      10 20 50
# 0  19 17 18
# 14 19 16 18
# 31  4 17  0
# 45  3  0  0
# 50  0 12  0       # unbalanced design

# type III SS
library(car)
ecklonia_anova <- aov(biomass ~ depth * time, data = ecklonia_decomp_21)
Anova(ecklonia_anova, type = "III")
#             Sum Sq  Df  F value    Pr(>F)    
# (Intercept) 19.9625   1 222.8055 < 2.2e-16 ***
# depth        0.1194   2   0.6662    0.5153    
# time         3.4022   1  37.9727 7.513e-09 ***
# depth:time   0.0655   2   0.3657    0.6944    
# Residuals   12.2747 137                       

ecklonia_20m_anova <- aov(biomass ~ site * time, data = ecklonia_mid)
Anova(ecklonia_20m_anova, type = "III")    
# Sum Sq Df  F value    Pr(>F)    
# (Intercept) 8.8264  1 103.3351 2.054e-14 ***
# site        0.0512  2   0.2997    0.7422    
# time        3.1138  1  36.4553 1.247e-07 ***
# site:time   0.0084  2   0.0489    0.9523    
# Residuals   4.8687 57 


ecklonia_14d <- subset(ecklonia_decomp, time=="14")
ecklonia_14d_anova <- aov(biomass ~ site, data = ecklonia_14d)
summary(ecklonia_14d_anova)



####### NESTED ANOVA ###############################################
# site needs to be nested in depth. set both as a facor
ecklonia_decomp_21$site <- as.factor(ecklonia_decomp_21$site)
ecklonia_nest_anova <- aov(biomass ~ depth/site, ecklonia_decomp_21, type="III")
summary(ecklonia_nest_anova)  
#               Df Sum Sq Mean Sq F value  Pr(>F)   
# depth         2  1.766  0.8828   5.394 0.00558 **
# epth:site    6  1.376  0.2293   1.401 0.21881   
# Residuals   134 21.931  0.1637 
TukeyHSD(ecklonia_nest_anova, "depth")
# $depth
# diff         lwr       upr     p adj
# 20-10 -0.08406511 -0.27183487 0.1037047 0.5398723
# 50-10  0.19358499 -0.02081308 0.4079831 0.0857199
# 50-20  0.27765009  0.07673997 0.4785602 0.0038192
table(ecklonia_decomp_21$site, ecklonia_decomp_21$depth)

####### BOX PLOT VISUALISATION ########################
library(minpack.lm)
library(ggpubr)
ggboxplot(ecklonia_decomp_21, x = "time", y = "biomass", color = "depth",
          palette = c("#00AFBB", "#E7B800", "green"))

########## EXPONENTIAL DECAY #########################################
# take the log(y) but not of the x's
ecklonia_decomp_summary_21_fit = lm(log(ecklonia_decomp_summary_21$mean) ~ ecklonia_decomp_summary_21$time)
summary(ecklonia_decomp_summary_21_fit)

factor=ecklonia_decomp_summary_21_fit$coefficients(ecklonia_decomp_summary_21$mean)
factor


########### GNLS MODEL #############################

# violin plot of gC remaining over time
decomp_2021 %>% 
  ggplot(aes(time, gC, group=time, fill=species)) +
  geom_violin() +
  facet_wrap(~species) +
  scale_fill_manual(values=c("#B09C85FF",
                             "#91D1C2FF")) +
  theme_bw()

decomp_2021_14 <- subset(decomp_2021, time == "14")

ggplot_boxplot(aes(depth, D, data=decomp_2021_14,fill=species))+
  facett_wrap(~species)+
  theme_bw()

decomp_2021 %>% 
  ggplot(aes(time, D, fill = depth, color = depth), na.rm = TRUE) +
  facet_wrap(~species) +
  geom_boxplot() +
  labs(x="Detrital age (days)", y= "% Carbon") +
  theme_bw() 



#### MORTEN EDITS MARCH 6 2023 #################

########## (%RC) - ECKLONIA SUMMARY DATA ########################
#################################################################

# ecklonia_decomp_21 <- subset (decomp_2021, species == "Ecklonia radiata")
# ecklonia_decomp_21$time <- factor(ecklonia_decomp_21$time)
# ecklonia_decomp_21$prC
# Run the functions mean, sd, and se on the value of "percent.biomass.remaining" for each group, 
# broken down by depth + time

ecklonia_decomp_summary_21 <- ddply (ecklonia_decomp_21, c("depth", "time"), summarise,
                                     N    = sum(!is.na(prC)),
                                     mean = mean(prC, na.rm=TRUE),
                                     sd   = sd(prC, na.rm = TRUE),
                                     se   = sd / sqrt(N))

ecklonia_decomp_summary_21
#   depth time  N       mean         sd         se
# 1    10    0  3 1.05306878 0.07802549 0.04504803
# 2    10   14 19 0.31155874 0.23131032 0.05306623
# 3    10   31  4 0.18281325 0.21373575 0.10686788
# 4    10   50  3 0.13574967 0.23512532 0.13574967
# 5    20   14 16 0.45421283 0.31697293 0.07924323
# 6    20   31 17 0.27760461 0.22718432 0.05510029
# 7    20   50 12 0.03515046 0.03437310 0.00992266
# 8    50   14 18 0.47893602 0.36052170 0.08497578
ecklonia_decomp_summary_21$time
#[1] 0  14 31 50 14 31 50 14

ecklonia_decomp_summary_21$mean
# [1] 1.05306878 0.31155874 0.18281325 0.13574967 0.45421283 0.27760461 0.03515046 0.47893602

# summarising data using summarySE fn w 95% ci
# need to rework: summarySE(ecklonia_decomp_21, measurevar="biomass", groupvars=c("depth","time"))
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

# DATA PREP
ecklonia_decomp_summary_21$depth <- as.factor(ecklonia_decomp_summary_21$depth)
# ecklonia_decomp_summary_21$time <- as.integer(ecklonia_decomp_summary_21$time)



########## (%RB) - SCYTOTHALIA SUMMARY DATA #####################
#################################################################

# scytothalia_decomp_21 <- subset (decomp_2021, species == "Scytothalia dorycarpa")
scyto_decomp_21$prC

# Run the functions mean, sd, and se on the value of "percent.biomass.remaining" for each group, 
# broken down by depth + time
# need to address NAs 

scyto_decomp_summary_21 <- ddply (scyto_decomp_21, c("depth", "time"), summarise,
                                        N    = sum(!is.na(prC)),
                                        mean = mean(prC, na.rm=TRUE),
                                        sd   = sd(prC, na.rm = TRUE),
                                        se   = sd / sqrt(N))
scyto_decomp_summary_21

# depth time  N      mean         sd         se
# 1    10    0 10 1.0000000 0.06306963 0.01994437
# 2    10   14 18 0.2294444 0.07981392 0.01881232
# 3    10   31  4 0.0975000 0.06238322 0.03119161
# 4    10   50  2 0.0950000 0.13435029 0.09500000
# 5    20   14 18 0.3422222 0.16611674 0.03915409
# 6    20   31 18 0.1588889 0.10317825 0.02431935
# 7    20   50 10 0.0710000 0.05404730 0.01709126
# 8    50   14 18 0.2305556 0.08495481 0.02002404

scytothalia_decomp_summary_21$time
#[1] 0  14 31 50 14 31 50 14

scytothalia_decomp_summary_21$mean
# [1] 1.0000000 0.2294444 0.0975000 0.0950000 0.3422222 0.1588889 0.0710000 0.23055561.0000000 0.2294444 0.0975000 0.0950000 0.3422222 0.1588889 0.0710000 0.2305556

# summarising data using summarySE fn w 95% ci
# need to rework: summarySE(scytothalia_decomp_21, measurevar="biomass", groupvars=c("depth","time"))
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

# DATA PREP
scyto_decomp_summary_21$depth <- factor(scyto_decomp_summary_21$depth)
# scytothalia_decomp_summary_21$time <- as.integer(scytothalia_decomp_summary_21$time)

###################################################################
# MODEL1EC - Ecklonia percent remainig carbon with depth ########
# START RUN HERE FOR L&O FIGURE 2B ######

decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$depth <- factor(decomp_2021$depth)
decomp_2021$species <- factor(decomp_2021$species)
decomp_2021$site <- factor(decomp_2021$site)
decomp_2021$gC <- (decomp_2021$C * decomp_2021$dw)/100
decomp_2021$D <- ((1-decomp_2021$prC)/decomp_2021$time)
ecklonia_decomp_21 <- subset(decomp_2021, species == "Ecklonia radiata")
#calculate mean initial gC for E.radiata
mean(1.4028406, 1.4260579, 1.6029657) #[1] 1.402841
ecklonia_decomp_21$prC <- (ecklonia_decomp_21$gC/1.402841)
ecklonia_decomp_21$prC[is.na(ecklonia_decomp_21$prC)] <- 0
ecklonia_decomp_21$prC #column created for percent remaining carbon

ecklonia_decomp_21$time #<- numeric(ecklonia_decomp_21$time)
table(ecklonia_decomp_21$time)
m1e <- nls(prC ~ SSasymp(time, A, I, logk), 
           data = ecklonia_decomp_21)
summary(m1e)
# Formula: prC ~ SSasymp(time, A, I, logk)
# Parameters:
#       Estimate Std. Error t value Pr(>|t|)    
# A     0.08953    0.07701   1.163    0.248    
# I     1.00897    0.15113   6.676 2.04e-09 ***
# logk -2.63187    0.28129  -9.356 6.87e-15 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2616 on 79 degrees of freedom

# NEW MODEL with values from m1 for start for each depth
# write out the entire eqn: simplify by removing I since it is 1
m2e <- nls(prC ~ exp(k[depth] * time),
           start = list(k = c(-exp(-2.63187),-exp(-2.63187),-exp(-2.63187))),
           data = ecklonia_decomp_21)
summary(m2e)
# Formula: prC ~ exp(k[depth] * time)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# k1 -0.079814   0.012544  -6.362 8.38e-09 ***
# k2 -0.049830   0.006249  -7.975 4.90e-12 ***
# k3 -0.052585   0.009127  -5.762 1.18e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2522 on 79 degrees of freedom

# m3e <- nls(prC ~ exp(k[depth] * time) + A[depth],
#            start = list(I = c(1.00897 ,1.00897 ,1.00897 ),
#                         k = c(-exp(-2.63187),-exp(-2.63187),-exp(-2.63187)),
#                         A = c(0.08953 ,0.08953 ,0.08953 )),
#            data = ecklonia_decomp_21)
# summary(m3e)

# ecklonia_decomp_21$time <- factor(ecklonia_decomp_21$time)
plot(resid(m2e)~ecklonia_decomp_21$time)
abline(0,0)
hist(resid(m2e))

#m3e model with Residual fraction (A) does not fit the model. Using gnls to normalize for time|depth
m4e <- gnls(prC ~ exp(k * time),
            start = list(k = c(-exp(-2.63187),-exp(-2.63187),-exp(-2.63187))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
plot(resid(m4e, type = "normalized") ~ time, data = ecklonia_decomp_21)
abline(0,0)
hist(resid(m4e, type = "normalized"))
summary(m4e) #ECKLONIA 10m (v 20m and 50m)
# Coefficients:
#   Value  Std.Error   t-value p-value
# k.(Intercept) -0.07608197 0.01323736 -5.747517  0.0000
# k.depth20      0.02576646 0.01400361  1.839987  0.0691
# k.depth50      0.02349709 0.01822170  1.289511  0.2006
intervals(m4e)
# Coefficients:
#   lower        est.       upper
# k.(Intercept) -0.10238432 -0.07608197 -0.04977961
# k.depth20     -0.00205842  0.02576646  0.05359135
# k.depth50     -0.01270904  0.02349709  0.05970322
# Variance function:
#   lower         est.        upper
# 10 -0.03106967 -0.016357950 -0.001646234
# 20 -0.03888385 -0.023937789 -0.008991724
# 50 -0.04034198 -0.004146673  0.032048636
# 
# Residual standard error:
#   lower      est.     upper 
# 0.2520553 0.3713059 0.5469754 

# change order of depth to get summary for 20m
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("20","50","10"))
m4e <- gnls(prC ~ exp(k * time),
            start = list(#I = c(1.00897, 1.00897, 1.00897),
                         k = c(exp(-2.63187),-exp(-2.63187),-exp(-2.63187))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
summary(m4e) #ECKLONIA 20m 
# Generalized nonlinear least squares fit
# Model: prC ~ exp(k * time) 
# Data: ecklonia_decomp_21 
# AIC      BIC   logLik
# 8.354131 26.00665 2.822934
# 
# Variance function:
#   Structure: Exponential of variance covariate, different strata
# Formula: ~time | depth 
# Parameter estimates:
#   10           20           50 
# -0.016357962 -0.023937793 -0.004146676 
# 
# Coefficients:
#   Value   Std.Error    t-value p-value
# k.(Intercept) -0.05031549 0.004568747 -11.012973  0.0000
# k.depth50     -0.00226939 0.013329521  -0.170253  0.8652
# k.depth10     -0.02576630 0.014003575  -1.839980  0.0691
# 
# Correlation: 
#   k.(In) k.dp50
# k.depth50 -0.343       
# k.depth10 -0.326  0.112
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.83100508 -0.70841015 -0.09637202  0.53988489  2.44752596 
# 
# Residual standard error: 0.377512 
# Degrees of freedom: 92 total; 89 residual
confint(m4e)
#                 2.5 %       97.5 %
# k.(Intercept) -0.05927007 -0.041360910
# k.depth50     -0.02839477  0.023855995
# k.depth10     -0.05321280  0.001680206
intervals(m4e)

# change order of depth to get summary for 50m
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("50","10","20"))
m4e <- gnls(prC ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
summary(m4e) #50m
# Generalized nonlinear least squares fit
# Model: prC ~ exp(k * time) 
# Data: ecklonia_decomp_21 
# AIC      BIC   logLik
# 8.354131 26.00665 2.822934
# 
# Variance function:
#   Structure: Exponential of variance covariate, different strata
# Formula: ~time | depth 
# Parameter estimates:
#   10           20           50 
# -0.016357747 -0.023937645 -0.004146369 
# 
# Coefficients:
#   Value  Std.Error   t-value p-value
# k.(Intercept) -0.05258488 0.01252209 -4.199370  0.0001
# k.depth10     -0.02349906 0.01822201 -1.289597  0.2005
# k.depth20      0.00226941 0.01332952  0.170255  0.8652
# 
# Correlation: 
#   k.(In) k.dp10
# k.depth10 -0.687       
# k.depth20 -0.939  0.646
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.83100977 -0.70840884 -0.09637292  0.53988499  2.44752425 
# 
# Residual standard error: 0.3775104 
# Degrees of freedom: 92 total; 89 residual
intervals(m4e)

###########################

ecklonia_decomp_21 <- subset(decomp_2021, species == "Ecklonia radiata")
#calculate mean initial gC for E.radiata
mean(1.4028406, 1.4260579, 1.6029657) #[1] 1.402841
ecklonia_decomp_21$prC <- (ecklonia_decomp_21$gC/1.402841)
ecklonia_decomp_21$prC[is.na(ecklonia_decomp_21$prC)] <- 0
ecklonia_decomp_21$prC #column created for percent remaining carbon
#ecklonia_decomp_21$time <- factor(ecklonia_decomp_21$time)
# m6e <- nls(prC ~ exp(k * time),
#            start = list(k = c(-exp(-2.63187)),
#                         data = ecklonia_decomp_21))
# summary (m6e)

ecklonia_new <- data.frame(time = rep(seq(0, 50, by = 1)),
                           depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new$depth <- factor(ecklonia_new$depth)
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth)
# depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new$fit <- predict(m4e, ecklonia_new)
#ecklonia_new$time <- factor(ecklonia_new$time)
predictions <- predict(m4e, ecklonia_new, interval = "confidence", level = 0.95)
predictions

ggplot() + geom_point(data = ecklonia_decomp_21, aes(time, prC, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=2), size =2, alpha = 0.3) +
  scale_color_manual(values=c("#00BFC4", "lightskyblue", "royalblue4")) +
  geom_line(data = ecklonia_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = ecklonia_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = ecklonia_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  labs(x="Days (in situ)", y="Proportion of initial carbon (dw)",  title= "Ecklonia radiata") +
  theme_bw()+
  ylim(0,1) +
  #geom_ribbon(data = predictions, aes(time, ymin = lwr, ymax = upr), alpha = 0.2) + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
ggsave("scyto_decomp_plot.pdf", width = 4, height = 4, units = "in")

# START RUN HERE FOR L&O EXPONENTIAL CARBON DECAY CALCULATION FOR EXPORT ESTIMATE ######
ecklonia_new_alldepths <- data.frame(time = rep(seq(0, 60, by = 1)),
                           depth = c(rep("10", 61), rep("20", 61), rep("50", 61)))
ecklonia_new_alldepths$fit <- predict(m4e, ecklonia_new_alldepths)
ecklonia_new_alldepths$depth <- factor(ecklonia_new_alldepths$depth)
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth)
install.packages("openxlsx")
library(openxlsx)
write.xlsx(ecklonia_new_alldepths, "Ecklonia POC export output.xlsx", rowNames = FALSE)
# depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new_alldepths$fit <- predict(m4e, ecklonia_new_alldepths)
ggplot() + geom_line(data = ecklonia_new_alldepths, aes(time, fit, color = depth))
###################################################################
###################################################################
# MODEL1S - scytothalia %rC with depth ################################

decomp_2021 <- read.csv ("2021 Marmion Kelp Litter Data - Raw.csv")
decomp_2021$depth <- factor(decomp_2021$depth)
decomp_2021$species <- factor(decomp_2021$species)
decomp_2021$site <- factor(decomp_2021$site)
decomp_2021$gC <- (decomp_2021$C * decomp_2021$dw)/100
decomp_2021$D <- ((1-decomp_2021$rb)/decomp_2021$time)
scyto_decomp_21 <- subset(decomp_2021, species == "Scytothalia dorycarpa")
scyto_decomp <- read.csv ("scyto_decomp.csv")
scyto_decomp$depth <- factor (scyto_decomp$depth)
scyto_decomp$time <- scyto_decomp$days
scyto_decomp$rb <- scyto_decomp$biomass
scyto_decomp_21$C
scyto_decomp_21$gC
#calculate mean initial gC for S.dorycarpa
mean(1.6968091, 1.9180787, 1.7150520, 1.5681666, 1.6064144, 1.6093168, 1.8729970, 1.6686021, 1.7789071, 1.6720632) #[1] 1.696809
scyto_decomp_21$prC <- (scyto_decomp_21$gC/1.696809)
scyto_decomp_21$prC[is.na(scyto_decomp_21$prC)] <- 0
scyto_decomp_21$prC #column created for percent remaining carbon
# scyto_decomp <- read.csv ("2021 Marmion Kelp Litter Data - Scyto.csv")
# scyto_decomp$species <- as.factor(scyto_decomp$species)
# scyto_decomp$depth <- as.factor(scyto_decomp$depth)
#scyto_decomp_21$time <- numeric(scyto_decomp_21$time)
summary(scyto_decomp_21)
scyto_decomp_21

m1s <- nls(prC ~ SSasymp(time, A, I, logk), data = scyto_decomp_21)
summary(m1s)
# Formula: rb ~ SSasymp(time, A, I, logk)
# Estimate Std. Error t value Pr(>|t|)    
# A     0.10504    0.02279   4.608 1.26e-05 ***
# I     1.00713    0.03467  29.048  < 2e-16 ***
# logk -2.09446    0.09491 -22.067  < 2e-16 ***

# for gnls
#       Value  Std.Error   t-value p-value
# A     0.1050392 0.02279293   4.608412       0
# I     1.0071298 0.03467157  29.047710       0
# logk -2.0944627 0.09491202 -22.067412       0

# NEW MODEL with values from m1 for start for each depth
# write out the entire eqn: simplify by removing I since it is 1
m2s <- nls(prC ~ exp(k[depth] * time),
           start = list(k = c(-exp(-2.09446 ),-exp(-2.09446 ),-exp(-2.09446))),
           data = scyto_decomp_21)
summary(m2s)
# Formula: biomass ~ exp(k[depth] * time)
# Parameters:
#     Estimate Std. Error t value Pr(>|t|)    
# k1 -0.104427   0.007567  -13.80   <2e-16 ***
# k2 -0.070157   0.003889  -18.04   <2e-16 ***
# k3 -0.105978   0.007856  -13.49   <2e-16 ***
plot(resid(m2s)~scyto_decomp_21$time)
abline(0,0)
hist(resid(m2s))

m3s <- nls(prC ~ (exp(k[depth] * time)) + A[depth],
           start = list(k = c(-exp(-2.09446 ),-exp(-2.09446 ),-exp(-2.09446)),
                        A = c(0.10504, 0.10504 ,0.10504)),
           data = scyto_decomp_21)
summary(m3s)
# Estimate Std. Error t value Pr(>|t|)    
# k1 -9.885e-02  1.025e-02  -9.644   <2e-16 ***
# k2 -6.503e-02  6.114e-03 -10.636   <2e-16 ***
# k3 -9.501e-02  9.914e-03  -9.583   <2e-16 ***
# A1  1.690e-02  2.492e-02   0.678    0.499    
# A2  1.675e-02  2.340e-02   0.716    0.476    
# A3  1.269e-05  2.711e-02   0.000    1.000   
plot(resid(m3s)~scyto_decomp$days)
abline(0,0)
hist(resid(m3s))


m4s <- gnls(rb ~ exp(k * time) + A,
            start = list(k = c(-exp(-2.1153834),-exp(-2.1153834),-exp(-2.1153834)),
                         A = c(0.1043025,0.1043025,0.1043025)),
            params = list(k ~ depth,
                          A ~ depth),
            weights = varExp(form = ~ time|depth),
            data = scyto_decomp)
plot(resid(m4s, type = "normalized") ~ time, data = scyto_decomp)
abline(0,0)
hist(resid(m4s, type = "normalized"))
summary(m4s) #SCYTOTHALIA 10m (v 20m and 50m)
# Coefficients:
# Value   Std.Error    t-value p-value
# k.(Intercept) -0.09739559 0.008545654 -11.397090  0.0000
# k.depth20      0.03281864 0.010203039   3.216555  0.0016
# k.depth50      0.00238330 0.012028928   0.198131  0.8433
# A.(Intercept)  0.01090224 0.020732640   0.525849  0.5999
# A.depth20     -0.00362188 0.029443493  -0.123011  0.9023
# A.depth50     -0.01088955 0.030119449  -0.361546  0.7183
intervals(m4s)


# calculate without A but with gnls
m4s <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = scyto_decomp_21)
summary(m4s) #k.(Intercept) for 10m data
# Value   Std.Error    t-value p-value
# k.(Intercept) -0.10270762 0.007650771 -13.424480  0.0000
# k.depth20      0.03410146 0.008649161   3.942748  0.0002
# k.depth50     -0.00209691 0.009799462  -0.213982  0.8310
# change order of depth to get summary for 20m

scyto_decomp_21$depth <- factor(scyto_decomp_21$depth, levels = c("20","50","10"))
m4s <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = scyto_decomp_21)
summary(m4s) #scyto 20m 
# Coefficients:
# Value   Std.Error    t-value p-value
# k.(Intercept) -0.06860776 0.004034210 -17.006491   0e+00
# k.depth50     -0.03619677 0.007332803  -4.936280   0e+00
# k.depth10     -0.03410018 0.008649182  -3.942591   2e-04
intervals(m4s)

# change order of depth to get summary for 50m
scyto_decomp_21$depth <- factor(scyto_decomp_21$depth, levels = c("50","10","20"))
m4s <- gnls(rb ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = scyto_decomp_21)
summary(m4s) #scyto 50m 
# Coefficients:
# k.(Intercept) -0.10480453 0.006123329 -17.115614   0.000
# k.depth10      0.00209688 0.009799458   0.213980   0.831
# k.depth20      0.03619824 0.007332731   4.936529   0.000
intervals(m4s)



# FIGURE 2.B ################################################

# reorder depth factor
scyto_decomp$depth <- factor(scyto_decomp$depth, levels = c("10","20","50"))

# Start by calculating means and SE for each depth prC (scyto)
scyto_decomp_summary_21 <- ddply (scyto_decomp_21, c("depth", "time"), summarise,
                               N    = sum(!is.na(rb)),
                               mean = mean(rb, na.rm=TRUE),
                               sd   = sd(rb, na.rm = TRUE),
                               se   = sd / sqrt(N))
scyto_decomp_summary_21
# make a fit line for decay estimated by model 4
scyto_new <- data.frame(time = rep(seq(0, 50, by = 1), 2),
                        depth = c(rep("10", 51), rep("20", 51)))
scyto_new$depth <- factor(scyto_new$depth)
scyto_new$fit <- predict(m4s, newdata =  scyto_new)
#decomp_new$time <- numeric(decomp_new$time)
scyto_new

# Plot prC means +-se for each species with decomp_new$fit
ggplot() + geom_point(data = scyto_decomp_21, aes(time, prC, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=3), size = 2, alpha = 0.2) +
  geom_line(data = scyto_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = scyto_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = scyto_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  scale_color_manual(values=c("#F8766D", "darkorange", "violet")) +
  labs(x="Days (in situ)", y="Proportion of initial Carbon mass (dw)",  title= "Scytothalia dorycarpa") +
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = 'bottom', legend.direction = "horizontal")


decomp_2021$time <- factor(decomp_2021$time)
model <- aov(prC~species*depth, data=decomp_2021)
summary(model)

decomp_2021$time <- factor(ecklonia_decomp_21$time)
model <- aov(prC~species*depth, data=ecklonia_decomp_21)
summary(model)

decomp_2021$time <- factor(decomp_2021$time)
model <- aov(prC~species*depth, data=decomp_2021)
summary(model)

# START RUN HERE FOR L&O EXPONENTIAL CARBON DECAY CALCULATION FOR EXPORT ESTIMATE -SCYTOTHALIA ######
scyto_new_alldepths <- data.frame(time = rep(seq(0, 60, by = 1)),
                                     depth = c(rep("10", 61), rep("20", 61), rep("50", 61)))
scyto_new_alldepths$fit <- predict(m4s, scyto_new_alldepths)
scyto_new_alldepths$depth <- factor(scyto_new_alldepths$depth)
scyto_decomp_21$depth <- factor(scyto_decomp_21$depth)
write.xlsx(scyto_new_alldepths, "scyto POC export output.xlsx", rowNames = FALSE)
# depth = c(rep("10", 51), rep("20", 51)))
scyto_new_alldepths$fit <- predict(m4s, scyto_new_alldepths)
ggplot() + geom_line(data = scyto_new_alldepths, aes(time, fit, color = depth))
























############# UNKN CPY PASTE ###### KEEP? #####
# change order of depth to get summary for 50m
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("50","10","20"))
m4e <- gnls(prC ~ exp(k * time),
            start = list(k = c(exp(-2.69769),-exp(-2.69769),-exp(-2.69769))),
            params = list(k ~ depth),
            weights = varExp(form = ~time|depth),
            data = ecklonia_decomp_21)
summary(m4e) #ECKLONIA 50m 
# Coefficients:
#   Value  Std.Error   t-value p-value
# k.(Intercept) -0.05317432 0.01193696 -4.454594  0.0000
# k.depth10     -0.02416154 0.01780970 -1.356650  0.1783
# k.depth20      0.00392236 0.01278317  0.306838  0.7597

# FIGURE 2.B ################################################

# reorder depth factor
ecklonia_decomp_21$depth <- factor(ecklonia_decomp_21$depth, levels = c("10","20","50"))

# Start by calculating means and SE for each depth prC (ECKLONIA)
ecklonia_decomp_summary_21 <- ddply (ecklonia_decomp_21, c("depth", "time"), summarise,
                                     N    = sum(!is.na(prC)),
                                     mean = mean(prC, na.rm=TRUE),
                                     sd   = sd(prC, na.rm = TRUE),
                                     se   = sd / sqrt(N))
ecklonia_decomp_summary_21
# make a fit line for decay estimated by model 4
ecklonia_new <- data.frame(time = rep(seq(0, 50, by = 1), 2),
                           depth = c(rep("10", 51), rep("20", 51)))
ecklonia_new$depth <- factor(ecklonia_new$depth)
ecklonia_new$fit <- predict(m4e, newdata =  ecklonia_new)
#decomp_new$time <- numeric(decomp_new$time)
ecklonia_new

# Plot prC means +-se for each species with decomp_new$fit
ggplot() + geom_point(data = ecklonia_decomp_21, aes(time, prC, color=depth, shape=depth), na.rm=TRUE, position=position_dodge(width=2), size =2, alpha = 0.2) +
  geom_line(data = ecklonia_new, aes(time, fit, color=depth), size=1, alpha=0.6) +
  geom_hline(yintercept=c(0.5), linetype="dotted", alpha = 0.5) +
  geom_errorbar(data = ecklonia_decomp_summary_21, aes(time, mean, ymin = mean - se, ymax = mean + se), na.rm = TRUE, width = 0.1) +
  geom_point(data = ecklonia_decomp_summary_21, aes(time, mean, color=depth, shape=depth), na.rm = TRUE, size = 3) + 
  scale_color_manual(values=c("#00BFC4", "lightskyblue", "royalblue4")) +
  labs(x="Days (in situ)", y="Proportion of initial biomass (dw)",  title= "Ecklonia radiata") +
  theme_bw()+
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.direction = "horizontal")









