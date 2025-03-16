#!/usr/bin/Rscript
# coding: utf-8

#***********************************************************************
# GBApy (Growth Balance Analysis for Python)
# Copyright © 2024-2025 Charles Rocabert
# Web: https://github.com/charlesrocabert/gbapy
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

library("tidyverse")
library("rstudioapi")
library("cowplot")
library("ggpmisc")
library("Matrix")
library("ggrepel")
library("latex2exp")


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)

d = read.table("../output/mmsyn_fcr_v1_state_optimum.csv", h=T, sep=";")

x      = c(0.0, 1e-05, 0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0)
mu_obs = c(0.42, 0.35, 0.32, 0.29, 0.3, 0.36, 0.3, 0.43, 0.3, 0.39, 0.38, 0.49, 0.47, 0.51, 0.51, 0.5, 0.52, 0.5, 0.46, 0.43)
mu_sim = d$mu
D      = data.frame(x, mu_obs, scale(mu_obs), mu_sim, scale(mu_sim))
names(D) = c("glc", "mu_obs", "mu_obs_scaled", "mu_sim", "mu_sim_scaled")

p1 = ggplot(D, aes(x, mu_obs)) +
  geom_point(aes(color="Observed")) +
  #geom_smooth(aes(color="Observed"), se=F, lty=2) +
  geom_line(data=D, aes(x, mu_sim, color="Simulated")) +
  scale_x_log10() +
  xlab("Added glucose (g/L)") +
  ylab("Growth rate") +
  ggtitle("Growth rate (absolute scale)") +
  theme_classic()
p2 = ggplot(D, aes(x, mu_obs_scaled)) +
  geom_point(aes(color="Observed")) +
  #geom_smooth(aes(color="Observed"), se=F, lty=2) +
  geom_line(data=D, aes(x, mu_sim_scaled, color="Simulated")) +
  scale_x_log10() +
  xlab("Added glucose (g/L)") +
  ylab("Re-scaled Growth rate") +
  ggtitle("Growth rate (re-scaled)") +
  theme_classic()
plot_grid(p1, p2, labels="AUTO")


