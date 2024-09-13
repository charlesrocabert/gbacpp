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

load_mass_fractions <- function()
{
  d = read.table("/Users/charlesrocabert/git/charlesrocabert/GBA_MMSYN/data/source/Breuer-et-al-2019/mass_fractions.csv", sep=";", h=T, check.names=F)
  rownames(d) = d$ID
  return(d)
}

load_proteomics <- function()
{
  d = read.table("/Users/charlesrocabert/git/charlesrocabert/GBA_MMSYN/data/fba_model/MMSYN_proteomics.csv", sep=";", h=T, check.names=F)
  return(d)
}

load_GPR <- function()
{
  d = read.table("/Users/charlesrocabert/git/charlesrocabert/GBA_MMSYN/data/fba_model/MMSYN_GPR.csv", sep=";", h=T, check.names=F)
  return(d)
}

load_trajectory_dmu <- function( conditions, confident )
{
  dmu = data.frame()
  for (condition in conditions)
  {
    df  = read.table(paste0("./output/MMSYN_",condition,"_state_trajectory.csv"), h=T, sep=";", check.names=F)
    dmu = rbind(dmu, df[dim(df)[1],])
  }
  dmu           = filter(dmu, condition%in%confident)
  dmu$condition = factor(dmu$condition, levels=unique(sort(as.integer(dmu$condition))))
  return(dmu)
}

load_trajectory_dp <- function( conditions, confident )
{
  dp = data.frame()
  for (condition in conditions)
  {
    df     = read.table(paste0("./output/MMSYN_",condition,"_p_trajectory.csv"), h=T, sep=";", check.names=F)
    df$P   = rowSums(df[,-which(names(df)%in%c("t","dt","index", "condition"))])
    df$phi = df$Ribosome/df$P
    dp     = rbind(dp, df[dim(df)[1],])
    
  }
  dp           = filter(dp, condition%in%confident)
  dp$condition = factor(dp$condition, levels=unique(sort(as.integer(dp$condition))))
  return(dp)
}

load_trajectory_db <- function( conditions, confident )
{
  db = data.frame()
  for (condition in conditions)
  {
    df = read.table(paste0("./output/MMSYN_",condition,"_b_trajectory.csv"), h=T, sep=";", check.names=F)
    db = rbind(db, df[dim(df)[1],])
    
  }
  db           = filter(db, condition%in%confident)
  db$condition = factor(db$condition, levels=unique(sort(as.integer(db$condition))))
  return(db)
}

build_mass_fraction_data <- function( db, mass_fractions, confident )
{
  M = data.frame()
  for(c in unique(db$condition))
  {
    df = filter(db, condition==c)
    df = df[dim(df)[1],]
    df = df[,-which(names(df)%in%c("t","dt","condition"))]
    df = df[,which(names(df)%in%mass_fractions$ID)]
    df = data.frame(names(df), t(df[1,]))
    names(df) = c("ID", "predicted")
    df$observed = mass_fractions[df$ID,"Fraction"]/100
    df$observed  = df$observed/sum(df$observed)
    df$condition = rep(c, dim(df)[1])
    df$rank      = rank(df$observed)
    M            = rbind(M, df)
  }
  M           = filter(M, condition%in%confident)
  M$condition = factor(M$condition, levels=unique(sort(as.integer(M$condition))))
  return(M)
}

build_proteomics_data <- function( dp, proteomics, gpr, confident )
{
  P = data.frame()
  for(c in unique(dp$condition))
  {
    ### Get last p line for condition c ###
    df = filter(dp, condition==c)
    df = df[dim(df)[1],]
    df = df[,-which(names(df)%in%c("t","dt","condition", "P", "phi"))]
    df  = df#/sum(df)
    ### Parse the GPR to build the proteome masses ###
    pf           = data.frame(proteomics$protein, rep(0.0, length(proteomics$protein)))
    names(pf)    = c("protein", "predicted")
    pf$observed  = proteomics$mass*0.001#/sum(proteomics$mass)
    rownames(pf) = pf$protein
    for(i in seq(1, dim(gpr)[1]))
    {
      reaction    = gpr$reaction[i]
      step1       = gpr$GPR[i]
      step2       = strsplit(step1,"|",fixed=T)
      total_stoic = 0.0
      #print(paste("> gpr =", step1, reaction))
      for(elmt in step2)
      {
        step3       = strsplit(elmt,":",fixed=T)
        protein     = step3[[1]][1]
        stoic       = as.numeric(step3[[1]][2])
        total_stoic = total_stoic+stoic
      }
      #print(paste("> stoic sum =", total_stoic))
      for(elmt in step2)
      {
        step3   = strsplit(elmt,":",fixed=T)
        protein = step3[[1]][1]
        stoic   = as.numeric(step3[[1]][2])
        #print(paste("> current stoic =", protein, stoic))
        if (reaction%in%names(df))
        {
          pf[protein,"predicted"] = pf[protein,"predicted"]+df[,reaction]*stoic/total_stoic
        }
      }
    }
    pf           = filter(pf, predicted > 0.0)
    pf           = filter(pf, observed > 0.0)
    pf$observed  = pf$observed/sum(pf$observed)
    pf$predicted = pf$predicted/sum(pf$predicted)
    pf$condition = rep(c, dim(pf)[1])
    P            = rbind(P, pf)
  }
  P             = filter(P, condition%in%confident)
  P$condition   = factor(P$condition, levels=unique(sort(as.integer(P$condition))))
  P$observedp3  = P$observed*3
  P$observedm3  = P$observed/3
  P$observedp10 = P$observed*10
  P$observedm10 = P$observed/10
  return(P)
}

##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/charlesrocabert/GBA_Evolution_2/")

#-------------------------------#
# 1) Load experimental datasets #
#-------------------------------#
confident      = seq(1,14)
#confident      = seq(1,30)
converged      = c(seq(1,14),seq(23,30))
GLC            = c(5.0, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0.078125, 0.0390625, 0.01953125, 0.009765625, 0.0048828125, 0.00244140625, 0.001220703125, 0.0006103515625, 0.00030517578125, 0.000152587890625, 7.62939453125e-05, 3.814697265625e-05, 1.9073486328125e-05, 9.5367431640625e-06, 4.76837158203125e-06, 2.384185791015625e-06, 1.1920928955078125e-06, 5.960464477539062e-07, 2.980232238769531e-07, 1.4901161193847656e-07, 7.450580596923828e-08, 3.725290298461914e-08, 1.862645149230957e-08, 9.313225746154785e-09)
mass_fractions = load_mass_fractions()
proteomics     = load_proteomics()
gpr            = load_GPR()

#-------------------------------#
# 2) Load simulation results    #
#-------------------------------#
dmu   = load_trajectory_dmu(c(seq(1,30)), confident)
dp    = load_trajectory_dp(c(seq(1,30)), confident)
db    = load_trajectory_db(c(seq(1,30)), confident)
D     = cbind(dmu, dp[,-which(names(dp)%in%c("t","dt","condition", "index"))])
D     = cbind(D, db[,-which(names(db)%in%c("t","dt","condition", "index"))])
D$glc = GLC[confident]

M = build_mass_fraction_data(db, mass_fractions, confident)
P = build_proteomics_data(dp, proteomics, gpr, confident)

# head(gpr)
# 
# names(D)
# 
# 
# ggplot(D, aes(glc, mu, color=condition)) +
#   scale_x_log10() +
#   geom_point() +
#   theme_classic() +
#   theme(legend.position="none")
# # 
# ggplot(D, aes(mu, phi)) +
#   geom_line() +
#   #scale_x_log10() +
#   #scale_y_log10() +
#   xlab("Growth rate") +
#   ylab("Phi (p_ribosome/sum(p))") +
#   theme_classic() +
#   theme(legend.position="none")

# ggplot(M, aes(observed, predicted)) +
#   facet_wrap(.~condition) +
#   geom_abline(slope=1, intercept=0, lty=2) +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("adj.R2", "p"))) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_classic() +
#   theme(legend.position="none")

ggplot(P, aes(log10(observed), log10(predicted))) +
  geom_line(data=P, aes(log10(observed), log10(observedp3)), col="black") +
  geom_line(data=P, aes(log10(observed), log10(observedm3)), col="black") +
  geom_line(data=P, aes(log10(observed), log10(observedp10)), col="black") +
  geom_line(data=P, aes(log10(observed), log10(observedm10)), col="black") +
  geom_abline(slope=1, intercept=0, lty=2) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2", "p"))) +
  geom_point() +
  facet_wrap(.~condition) +
  theme_classic()
