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
library("renz")

### Load mass fractions data ###
load_mass_fractions <- function()
{
  d           = read.table("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/manual_curation/MMSYN_mass_fractions.csv", sep=";", h=T, check.names=F)
  rownames(d) = d$ID
  return(d)
}

### Build the mass fraction data ###
build_mass_fractions_data <- function( d_b, i )
{
  d_mf     = load_mass_fractions()
  D        = d_b[i,]
  D        = D[,-which(names(D)%in%c("condition", "iter", "t", "dt", "h2o"))]
  D        = D[,which(names(D)%in%d_mf$ID)]
  D        = data.frame(names(D), t(D))
  names(D) = c("id", "sim")
  D$obs    = d_mf[D$id,"Fraction"]
  D$obs    = D$obs/sum(D$obs)
  D$sim    = D$sim/sum(D$sim)
  D$obsp3  = D$obs*3
  D$obsm3  = D$obs/3
  D$obsp10 = D$obs*10
  D$obsm10 = D$obs/10
  return(D)
}

### Correlation of mass fraction prediction ###
mass_fractions_cor <- function( d_b, i )
{
  D     = build_mass_fractions_data(d_b, i)
  reg   = lm(log10(D$sim)~log10(D$obs))
  r2    = summary(reg)$adj.r.squared
  pval  = summary(reg)$coefficients[,4][[2]]
  SSres = sum((log10(D$sim)-log10(D$obs))^2)
  M     = mean(log10(D$obs))
  SStot = sum((log10(D$obs)-M)^2)
  R2    = 1-SSres/SStot
  return(c(r2, pval, R2))
}

### Load proteomics data ###
load_proteomics <- function()
{
  d           = read.table(paste0("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/manual_curation/MMSYN_proteomics.csv"), sep=";", h=T, check.names=F)
  d$protein   = str_replace(d$locus, "JCVISYN3A", "protein")
  d$mass      = d$mg_per_gP*0.001
  rownames(d) = d$protein
  d           = filter(d, mass > 0.0)
  return(d)
}

### Load protein contributions ###
load_protein_contributions <- function( model_path, model_name )
{
  filename = paste0(model_path,"/",model_name,"/protein_contributions.csv")
  d        = read.table(filename, h=T, sep=";", check.names=F)
  return(d)
}

### Load enzyme composition ###
# Note:
# Some reactions are filtered out as they relate to
# the model size reduction.
load_enzyme_composition <- function()
{
  AMINO_ACID_TRANSPORTERS = c("ALAt2r", "ARGt2r", "ASNt2r", "ASPt2pr", "CYSt2r", "GLNt2r",
                              "GLUt2pr", "GLYt2r", "HISt2r", "ISOt2r", "LEUt2r", "LYSt2r",
                              "METt2r", "PHEt2r", "PROt2r", "SERt2r", "THRt2r", "TRPt2r",
                              "TYRt2r", "VALt2r")
  CHARGING_REACTIONS = c("ALATRS", "ARGTRS", "ASNTRS", "ASPTRS", "CYSTRS", "GLUTRS",
                         "GLYTRS", "HISTRS", "ILETRS", "LEUTRS", "LYSTRS", "METTRS",
                         "PHETRS", "PROTRS", "SERTRS", "THRTRS", "TRPTRS", "TYRTRS",
                         "VALTRS")
  filename = paste0("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/manual_curation/MMSYN_enzyme_composition.csv")
  d        = read.table(filename, h=T, sep=";", check.names=F)
  d        = filter(d, !info_sample_rid %in% c("Protein_transl", AMINO_ACID_TRANSPORTERS, CHARGING_REACTIONS))
  return(d)
}

### Calculate simulated proteomics ###
calculate_simulated_proteomics <- function( d_p, protein_contributions, i )
{
  X             = d_p[i, -which(names(d_p)%in%c("condition", "iter", "t", "dt"))]
  r_ids         = names(X)
  p_ids         = unique(protein_contributions$protein)
  res           = data.frame(p_ids, rep(0.0, length(p_ids)))
  names(res)    = c("p_id", "value")
  rownames(res) = res$p_id
  for (r_id in r_ids)
  {
    e_conc  = X[r_id][[1]]
    contrib = protein_contributions[protein_contributions$reaction==r_id,]
    if (dim(contrib)[1] > 0)
    {
      for (j in seq(1, dim(contrib)[1]))
      {
        p_id              = contrib$protein[j]
        value             = e_conc*contrib$contribution[j]
        res[p_id,"value"] = res[p_id,"value"]+value
      }
    }
  }
  res$value = res$value/sum(res$value)
  return(res)
}

### Collect proteins excluded from the model ###
collect_excluded_proteins <- function( list_of_reactions, enzyme_composition )
{
  excluded_proteins = c()
  for(i in seq(1, dim(enzyme_composition)[1]))
  {
    r_id = enzyme_composition$info_sample_rid[i]
    if (!r_id %in% list_of_reactions | r_id %in% c("Ribosome", "AAabc", "AATRS"))
    {
      step1 = enzyme_composition$composition[i]
      step2 = strsplit(step1,"|",fixed=T)
      for (elmt in step2[[1]])
      {
        step3             = strsplit(elmt,",",fixed=T)
        prot_id           = str_replace(step3[[1]][1], "gene=JCVISYN3A_", "")
        prot_id           = str_replace(prot_id, " ", "")
        excluded_proteins = c(excluded_proteins, paste0("protein_",prot_id))
      }
    }
  }
  return(excluded_proteins)
}

### Build the proteome fraction data ###
build_proteomics_data <- function( d_p, i )
{
  obs_proteomics = load_proteomics()
  pcontrib       = load_protein_contributions(model_path, model_name)
  enz_comp       = load_enzyme_composition()
  reaction_ids   = names(d_p[,-which(names(d_p)%in%c("condition", "iter", "t", "dt"))])
  excluded       = collect_excluded_proteins(reaction_ids, enz_comp)
  sim_proteomics = calculate_simulated_proteomics(d_p, pcontrib, i)
  D              = data.frame(sim_proteomics$p_id, sim_proteomics$value)
  names(D)       = c("id", "sim")
  obs_sum        = sum(obs_proteomics$mass)
  sim_sum        = sum(D$sim)
  D              = filter(D, id%in%obs_proteomics$protein)
  D$obs          = obs_proteomics[D$id, "mass"]
  D              = D[!D$id%in%excluded,]
  D$obs          = D$obs/sum(D$obs)
  D$sim          = D$sim/sum(D$sim)#*0.23#/sim_sum#*0.23
  D$obsp3        = D$obs*3
  D$obsm3        = D$obs/3
  D$obsp10       = D$obs*10
  D$obsm10       = D$obs/10
  
  return(D)
}

### Correlation of proteomics prediction ###
proteomics_cor <- function( d_p )
{
  D     = build_proteomics_data(d_p, dim(d_p)[1])
  reg   = lm(log10(D$sim)~log10(D$obs))
  r2    = summary(reg)$adj.r.squared
  pval  = summary(reg)$coefficients[,4][[2]]
  SSres = sum((log10(D$sim)-log10(D$obs))^2)
  M     = mean(log10(D$obs))
  SStot = sum((log10(D$obs)-M)^2)
  R2    = 1-SSres/SStot
  return(c(r2, pval, R2))
}


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)

### Build the glucose vector ###
x        = c(0.0, 1e-05, 0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0)
x_adjust = x+0.4
x_raw    = x[2:length(x)]
x_raw    = x_raw[-3]
glc      = c(x_adjust, x_raw)

### Build the experimental dataset ###
mu_obs          = c(0.42, 0.35, 0.32, 0.29, 0.3, 0.36, 0.3, 0.43, 0.3, 0.39, 0.38, 0.49, 0.47, 0.51, 0.51, 0.5, 0.52, 0.5, 0.46, 0.43)
cond            = seq(1,20)
d_obs           = data.frame(cond, x_adjust, mu_obs)
names(d_obs)    = c("condition", "glc", "mu")

d_state = read.table("../output/mmsyn_fcr_v1_state_optimum.csv", h=T, sep=";")
d_b     = read.table("../output/mmsyn_fcr_v1_b_optimum.csv", h=T, sep=";")

d_state$glc    = glc
d_state$mu_obs = d_obs[d_state$condition, "mu"]
d_b$glc        = glc

MF_cor = data.frame()
for (i in seq(1, dim(d_b)[1]))
{
  res     = mass_fractions_cor(d_b, i )
  MF_cor = rbind(MF_cor, res)
}
MF_cor$condition = d_b$condition
MF_cor$glc       = glc
names(MF_cor) = c("r2", "pval", "R2", "condition", "glc")


ggplot(MF_cor, aes(glc, R2)) +
  geom_line() +
  scale_x_log10()

ggplot(d_state, aes(glc, mu)) +
  geom_line() +
  geom_point(aes(glc, mu_obs, color="Observed")) +
  scale_x_log10() +
  #scale_y_log10() +
  theme_classic()




ggplot(d, aes(glc, mu)) +
  geom_hline(yintercept=0.34) +
  geom_line() +
  geom_point(aes(glc, mu_obs, color="Observed")) +
  scale_x_log10() +
  #scale_y_log10() +
  xlab("Added glucose (g/L)") +
  ylab("Growth rate") +
  ggtitle("Growth rate (absolute scale)") +
  theme_classic()

# d_obs$log_glc = log10(d_obs$glc)
# dir.MM(d_obs[d_obs$glc<6, c("glc", "mu")], unit_v="Growth rate")
# d_obs
# res = dir.MM(d_obs[d_obs$glc<2, c("glc", "mu")], unit_v="Growth rate")
# 
# plot(res$data$S, res$data$v, pch=20, xlim=c(0,1.4), ylim=c(0, 0.6))
# lines(res$data$S, res$data$fitted_v)

