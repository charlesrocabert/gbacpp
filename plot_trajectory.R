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

##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)

d1 = read.table("./output/MMSYN_state_trajectory.csv", h=T, sep=";")
d2 = read.table("./output/MMSYN_b_trajectory.csv", h=T, sep=";")
d3 = read.table("./output/MMSYN_p_trajectory.csv", h=T, sep=";")
d4 = read.table("./output/MMSYN_c_trajectory.csv", h=T, sep=";")
m  = read.table("../GBA_MMSYN/data/source/Breuer-et-al-2019/MMSYN_mass_fractions.csv", h=T, sep=";", check.names=F)

d3$P   = rowSums(d3[,-which(names(d3)%in%c("t","dt","index"))])
d3$phi = d3$Ribosome/d3$P
d3$mu  = d1$mu

d1$index = seq(1, dim(d1)[1])*100
d2$index = seq(1, dim(d2)[1])*100
d3$index = seq(1, dim(d3)[1])*100

d1$t[1] = d1$t[2]
reg = lm(log10(d1$mu)~log10(d1$t))
d1$fit = 10^reg$coefficients[[1]]*d1$t^reg$coefficients[[2]]

p1 = ggplot(d1, aes(t, mu)) +
  geom_line(aes(t, fit), col="red", lty=2) +
  geom_line() +
  #scale_y_log10() +
  theme_classic() 
p2 = ggplot(d1, aes(index, log10(dt))) +
  geom_smooth(se=F) +
  geom_line(lty=2) +
  theme_classic()
p3 = ggplot(d1, aes(index, log10(mu_diff))) +
  geom_smooth(se=F) +
  geom_hline(yintercept=-10, color="red") +
  geom_line(lty=2) +
  theme_classic()
p4 = ggplot(d2, aes(t, Protein*100)) +
  geom_line() +
  theme_classic()
p5 = ggplot(d3, aes(mu, phi)) +
  geom_line() +
  theme_classic()

plot_grid(p1, p2, p3, p4, p5, ncol=2)


# df2 = d2 %>% rowwise() %>% pivot_longer(-t)
# 
# df2 = filter(df2, !name%in%c("dt", "index"))
# p4 = ggplot(df2, aes(t, value, color=name)) +
#   geom_line() +
#   ggtitle("Mass fractions") +
#   theme_classic() +
#   theme(legend.position="none")
# p4

dim(d1)
