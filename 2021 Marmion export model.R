# Particle tracking model estimates for particle export (days) from source to deep ocean (200m)

setwd("~/Desktop/R_litterbag")
particle_MAR <- read.csv("export estimates 0-60 days.csv")
particle_MAR$depth <- as.factor(particle_MAR$source.reef.depth)
particle_MAR$time <- particle_MAR$Particle.age..days.
particle_MAR$export <- particle_MAR$percentage.of.particles.past.shelf.break..200m.

library(tidyverse)
library(ggplot2)


p <- particle_MAR %>% 
  #na.omit() %>% 
  ggplot(aes(time, export,
             color = depth))+
  scale_color_manual(values=c("blue", "darkgreen")) +
  #geom_point(size = 2.5, alpha = 0.5) + # alpha is transparency +
  geom_line(size=1.5) +
  #geom_smooth() + # if you want a linear model, say method = lm in parenthesese
  geom_vline(xintercept=c(14,30), linetype="dotted", alpha = 0.5) +
  #geom_text(aes(x=14, label="2 weeks", y=92), colour="black", angle=90, vjust = 0, text=element_text(size=11)) +
  #geom_text(aes(x=30, label="1 month", y=92), colour="black", angle=90, vjust = 0) +
  ylim(0,100) +
  scale_x_continuous(breaks=seq(0,100,10)) + 
  scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "% of particles exported offshore (> 200 m)")) +
  guides(y = "none")+
  labs(title = "% of kelp particulate exported to deep ocean (source reef to shelf-break (200m))")+
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(text = element_text(family = "Times New Roman"))
p
p_wrapped <- wrap_labs(p, "normal")

p + g
#WITH PRED CURVE OF POC DECOMP
ggplot() + geom_point(data = particle_MAR, aes(time, export, color = depth), size = 2.5, alpha = 0.5) +
  scale_color_manual(values=c("blue", "darkgreen")) +
  geom_line() +
  geom_line(data = decomp_new, aes(time, fit, color=species), size=1, alpha=0.6) +
  scale_y_continuous(
    labels = scales::comma,
    name = "% remaining POC",
    sec.axis = sec_axis(decomp_new$time, name="% particles to deep (200m)", labels = scales::comma)) +
  #geom_smooth() + # if you want a linear model, say method = lm in parenthesese
  geom_vline(xintercept=c(14,30), linetype="dotted", alpha = 0.5) +
  #geom_text(aes(x=14, label="2 weeks", y=92), colour="black", angle=90, vjust = 0, text=element_text(size=11)) +
  #geom_text(aes(x=30, label="1 month", y=92), colour="black", angle=90, vjust = 0) +
  ylim(0,100) +
  scale_x_continuous(breaks=seq(0,100,10)) + 
  labs(title = "% of kelp particulate exported to deep ocean (source reef to shelf-break (200m))")+
  theme_bw()

# Particle tracking model estimates for particle export (days) from source to deep ocean (200m)

poc_MAR <- read.csv("POC Export Calculations - CSV format.csv")
poc_MAR$depth <- as.factor(poc_MAR$depth)
poc_MAR$time <- poc_MAR$time
poc_MAR$export <- poc_MAR$poc_exp
# Subset export data for Ecklonia radiata
poc_MAR$export_ecklonia <- poc_MAR$export[poc_MAR$species == "Ecklonia radiata"]
# Subset export data for Scytothalia dorycarpa
poc_MAR$export_scyto <- poc_MAR$export[poc_MAR$species == "Scytothalia dorycarpa"]
# Calculate the mean of export_ecklonia and export_scyto
poc_MAR$export_allsp <- rowMeans(poc_MAR[, c("export_ecklonia", "export_scyto")], na.rm = TRUE)
# Print the first few rows of the updated data frame
head(poc_MAR)

library(tidyverse)
library(ggplot2)

p <- poc_MAR %>% 
  ggplot(aes(time, color = depth)) +
  geom_ribbon(aes(ymin = export_ecklonia, ymax = export_scyto, fill = factor(depth)), alpha = 0.5, color = NA) +
  geom_vline(xintercept = 20, linetype = "dotted", alpha = 0.5) +
  ylim(0, 100) +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 30)) + 
  scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "% of POC exported offshore (> 200 m)")) +
  scale_color_manual(values = c("50" = "darkgreen", "15" = "royalblue")) +
  scale_fill_manual(values = c("50" = "darkgreen", "15" = "royalblue"), 
                    breaks = c(15, 50), labels = c("Depth 15", "Depth 50")) +
  guides(y = "none", color = guide_legend(title = "Depth")) +
  labs(title = "% of kelp particulate exported to deep ocean (source reef to shelf-break (200m))") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal",
        text = element_text(family = "Times New Roman"))
p

p <- poc_MAR %>% 
  ggplot(aes(time, color = depth)) +
  geom_ribbon(aes(ymin = export_ecklonia, ymax = export_scyto, fill = factor(depth)), alpha = 0.5) +
  geom_vline(xintercept = 20, linetype = "dotted", alpha = 0.5) +
  ylim(0, 100) +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 30)) + 
  scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "% of POC exported offshore (> 200 m)")) +
  scale_color_manual(values = c("50" = "darkgreen", "15" = "royalblue")) +
  scale_fill_manual(values = c("50" = "darkgreen", "15" = "royalblue"), 
                    breaks = c(15, 50), labels = c("Depth 15", "Depth 50")) +
  guides(y = "none", color = guide_legend(title = "Depth")) +
  labs(title = "% of kelp particulate exported to deep ocean (source reef to shelf-break (200m))") +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.direction = "horizontal",
        text = element_text(family = "Times New Roman"))
p


p + g
#WITH PRED CURVE OF POC DECOMP
ggplot() + geom_point(data = poc_MAR, aes(time, export, color = depth), size = 2.5, alpha = 0.5) +
  scale_color_manual(values=c("blue", "darkgreen")) +
  geom_line() +
  geom_line(data = decomp_new, aes(time, fit, color=species), size=1, alpha=0.6) +
  scale_y_continuous(
    labels = scales::comma,
    name = "% remaining POC",
    sec.axis = sec_axis(decomp_new$time, name="% pocs to deep (200m)", labels = scales::comma)) +
  #geom_smooth() + # if you want a linear model, say method = lm in parenthesese
  geom_vline(xintercept=c(14,30), linetype="dotted", alpha = 0.5) +
  #geom_text(aes(x=14, label="2 weeks", y=92), colour="black", angle=90, vjust = 0, text=element_text(size=11)) +
  #geom_text(aes(x=30, label="1 month", y=92), colour="black", angle=90, vjust = 0) +
  ylim(0,100) +
  scale_x_continuous(breaks=seq(0,100,10)) + 
  labs(title = "% of kelp particulate exported to deep ocean (source reef to shelf-break (200m))")+
  theme_bw()

