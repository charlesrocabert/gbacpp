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

d1 = read.table("./output/EC12b_model_optimum.csv", h=T, sep=";")
d2 = read.table("./output/EC12b_f_optimum.csv", h=T, sep=";")

p1 = ggplot(d, aes(condition, mu)) +
  geom_line() +
  ggtitle("Growth rate") +
  theme_classic()

df2 = d2 %>% rowwise() %>% pivot_longer(-condition)
p2 = ggplot(df2, aes(condition, value, color=name)) +
  geom_line() +
  ggtitle("Fluxes") +
  theme_classic()

plot_grid(p1, p2, ncol=1)

