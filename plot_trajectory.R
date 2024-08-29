#!/usr/bin/Rscript
# coding: utf-8

#***************************************************************************
# Copyright © 2024 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBA_MMSYN
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

d = read.table("./output/model_trajectory.csv", h=T, sep=";")

p1 = ggplot(d, aes(t, mu)) + geom_line() + theme_classic()

