#!/usr/bin/Rscript
# coding: utf-8

#***************************************************************************
# Copyright © 2024 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBA_Evolution_2
#
# plot_trajectory.py
# ------------------
# Plot the current optimization trajectory
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("rstudioapi")
library("cowplot")
library("ggpmisc")

compare_mass_fractions <- function( observed_m, c_trajectory, index )
{
  observed_m = filter(observed_m, ID%in%names(c_trajectory))
  observed_m$Fraction = observed_m$Fraction/sum(observed_m$Fraction)*100
  #print(observed_m)
  #stop()
  X        = c_trajectory[index,-which(names(c_trajectory)%in%c("condition","t","dt", "h2o"))]
  H2O      = c_trajectory[index,"h2o"]
  Xsum     = sum(t(X))
  X        = X/Xsum
  X        = t(X[observed_m$ID])
  D        = data.frame(observed_m$ID, observed_m$Fraction, X*100)
  names(D) = c("ID", "observed", "predicted")
  D        = D[order(D$observed, decreasing=T),]
  D$ID     = factor(D$ID, levels=D$ID)
  D$index  = seq(1, dim(D)[1])
  return(D)
}

##################
#      MAIN      #
##################

#directory = dirname(getActiveDocumentContext()$path)
#setwd(directory)
setwd("/Users/charlesrocabert/git/charlesrocabert/GBA_Evolution_2/")


path       = "/Users/charlesrocabert/Desktop/mmsyn_output"
model_name = "mmsyn_test"
condition  = "1"

d1 = read.table(paste0(path,"/",model_name,"_",condition,"_state_trajectory.csv"), h=T, sep=";")
d2 = read.table(paste0(path,"/",model_name,"_",condition,"_b_trajectory.csv"), h=T, sep=";")
d3 = read.table(paste0(path,"/",model_name,"_",condition,"_p_trajectory.csv"), h=T, sep=";")
d4 = read.table(paste0(path,"/",model_name,"_",condition,"_f_trajectory.csv"), h=T, sep=";")
d5 = read.table(paste0(path,"/",model_name,"_",condition,"_c_trajectory.csv"), h=T, sep=";")
m  = read.table("../GBA_MMSYN/data/source/Breuer-et-al-2019/mass_fractions.csv", h=T, sep=";", check.names=F)
p  = read.table("../GBA_MMSYN/data/source/Peter-MMSYN/MMSYN_proteomics.csv", h=T, sep=";", check.names=F)
names(d2)
X = t(d4[dim(d4)[1],])
names(X) = names(d4)
X = X[3:length(X)]
X = sort(X)

#barplot(log10(X), las=2, cex.names=0.5)

d3$P   = rowSums(d3[,-which(names(d3)%in%c("t","dt","index"))])
d3$phi = d3$Ribosome/d3$P
d3$mu  = d1$mu

d1$index = seq(1, dim(d1)[1])*500
d2$index = seq(1, dim(d2)[1])*500
d3$index = seq(1, dim(d3)[1])*500

d1$t[1] = d1$t[2]
reg = lm(log10(d1$mu)~log10(d1$t))
d1$fit = 10^reg$coefficients[[1]]*d1$t^reg$coefficients[[2]]

p1 = ggplot(d1, aes(index, mu)) +
  geom_line(aes(index, fit), col="red", lty=2) +
  geom_line() +
  #scale_y_log10() +
  theme_classic() 
p2 = ggplot(d1, aes(index, log10(dt))) +
  geom_smooth(se=F) +
  geom_line(lty=2) +
  theme_classic()
p3 = ggplot(d1, aes(index, log10(mu_diff))) +
  geom_smooth(se=F) +
  geom_hline(yintercept=-10, color="red", lty=2) +
  geom_line(lty=2) +
  theme_classic()
p4 = ggplot(d2, aes(t, Protein*100)) +
  geom_line() +
  theme_classic()
p5 = ggplot(d3, aes(mu, phi)) +
  geom_line() +
  theme_classic()

D = compare_mass_fractions(m, d5, dim(d5)[1])
D["Protein",]
p6 = ggplot(D, aes(ID, log10(predicted))) +
  geom_bar(stat="identity") +
  geom_point(aes(ID, log10(observed)), col="red") +
  #geom_smooth(se=F, method="lm") +
  #geom_smooth(aes(index, log10(observed)), se=F, method="lm") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p7 = ggplot(D, aes(observed, predicted)) +
  geom_smooth(method="lm", se=F) +
  geom_abline(slope=1, intercept=0, lty=2) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic()

mf_lss  = c()
mf_r2   = c()
mf_pval = c()
for(i in seq(1, dim(d2)[1]))
{
  df = compare_mass_fractions(m, d2, i)
  df$rank = seq(1, dim(df)[1])
  lss = sum((log10(df$observed)-log10(df$predicted))^2)/dim(df)[1]
  #reg = lm(log10(df$observed)~df$rank)
  reg = lm(log10(df$observed)~log10(df$predicted))
  r2 = summary(reg)$adj.r.squared
  pval = summary(reg)$coefficients[2,4]
  mf_lss = c(mf_lss, lss)
  mf_r2 = c(mf_r2, r2)
  mf_pval = c(mf_pval, pval)
}
d1$mf_lss = mf_lss
d1$mf_r2  = mf_r2
d1$mf_pval  = mf_pval

p8 = ggplot(d1, aes(index, log10(mf_lss))) +
  geom_line() +
  theme_classic()

p9 = ggplot(d1, aes(index, (mf_r2))) +
  geom_line() +
  theme_classic()

p10 = ggplot(d1, aes(index, log10(mf_pval))) +
  geom_line() +
  geom_hline(yintercept=log10(0.05), col="red", lty=2) +
  theme_classic()

p = plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=2)
ggsave(paste0("trajectory_",condition,".pdf"), p, width=15/2, height=15/2)

p
# df2 = d2 %>% rowwise() %>% pivot_longer(-t)
# 
# df2 = filter(df2, !name%in%c("dt", "index"))
# p4 = ggplot(df2, aes(t, value, color=name)) +
#   geom_line() +
#   ggtitle("Mass fractions") +
#   theme_classic() +
#   theme(legend.position="none")
# p4

# par(mfrow=c(1,2))
# plot(d2$DNA, type="l")
# plot(d2$RNA, type="l")

# d = read.table(paste0(path,"/",model_name,"_",condition,"_c_trajectory.csv"), h=T, sep=";")
# par(mfrow=c(3,1))
# plot(d$h2o, type="l")
# plot(d$DNA, type="l")
# plot(d$clpn, type="l")
# 
# d = d[,-which(names(d)%in%c("condition","t","dt", "h2o"))]
# d$B = rowSums(d)
# plot(d$Protein/d$B)
