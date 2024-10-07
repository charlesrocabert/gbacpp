#!/usr/bin/Rscript
# coding: utf-8

#***************************************************************************
# Copyright © 2024 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBA_Evolution_2
#
# optimum_results.R
# ------------------
# Plot optimum analysis results.
# (LOCAL SCRIPT)
#***************************************************************************

library("rstudioapi")
library("tidyverse")
library("cowplot")
library("RColorBrewer")
library("plotly")

##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/Desktop/fcr_output")

#-------------------------------#
# 1) Load experimental datasets #
#-------------------------------#
OPT = read.table("FCR_EFM2_state_optimum.csv", sep=";", h=T)
opt2 = read.table("FCR_EFM2_f_optimum.csv", sep=";", h=T)
OPT  = cbind(OPT, opt2[,-which(names(opt2)%in%c("condition", "t", "dt"))])
D = data.frame()
for(i in seq(1,100))
{
  d1 = read.table(paste0("FCR_EFM2_",i,"_state_trajectory.csv"), sep=";", h=T)
  d2 = read.table(paste0("FCR_EFM2_",i,"_f_trajectory.csv"), sep=";", h=T)
  d  = cbind(d1, d2[,-which(names(d2)%in%c("condition", "t", "dt"))])
  D = rbind(D, d)
}
p1 = ggplot(D, aes(rxn3, rxn4)) +
  geom_line(color="lightgrey") +
  geom_point(data=OPT, aes(rxn3, rxn4), color="red", size=3) +
  #geom_hline(yintercept=0.984174/2) +
  #geom_vline(xintercept=0.984174/2) +
  #geom_abline(intercept=0.984174, slope=-1) +
  #ylim(0.6, max(D$mu)) +
  theme_classic() +
  theme(legend.position="none")
p2 = ggplot(D, aes(rxn3, mu)) +
  geom_line() +
  geom_point(data=OPT, aes(rxn3, mu), color="red", size=3) +
  ylim(0.5, max(D$mu)) +
  theme_classic() +
  theme(legend.position="none")
plot_grid(p1, p2, ncol=1)

#p = plot_ly(x=D$rxn3, y=D$rxn4, z=D$mu, type="scatter3d", mode="lines", color=factor(D$condition))
#add_trace(p, x=OPT$rxn3, y=OPT$rxn4, z=OPT$mu, type="scatter3d", mode="markers", color=OPT$condition)
