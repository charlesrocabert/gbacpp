#!/usr/bin/Rscript
# coding: utf-8

#***************************************************************************
# Copyright © 2024 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBA_Evolution_2
#
# MMSYN_results.R
# ---------------
# Plot MMSYN results.
# (LOCAL SCRIPT)
#***************************************************************************

library("rstudioapi")
library("tidyverse")
library("cowplot")
library("RColorBrewer")

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

build_mass_fraction_data <- function( db, mass_fractions )
{
  M = data.frame()
  for(c in unique(db$condition))
  {
    df = filter(db, condition==c)
    df = df[dim(df)[1],]
    df = df[,-which(names(df)%in%c("condition", "glc"))]
    df = df[,which(names(df)%in%mass_fractions$ID)]
    df = data.frame(names(df), t(df[1,]))
    names(df)    = c("ID", "predicted")
    df$observed  = mass_fractions[df$ID,"Fraction"]/100
    #df$predicted = df$predicted/sum(df$predicted)
    #df$observed  = df$observed/sum(df$observed)
    df$condition = rep(c, dim(df)[1])
    df$rank      = rank(df$observed)
    M            = rbind(M, df)
  }
  M$condition   = factor(M$condition, levels=unique(sort(as.integer(M$condition))))
  M$observedp3  = M$observed*3
  M$observedm3  = M$observed/3
  M$observedp10 = M$observed*10
  M$observedm10 = M$observed/10
  return(M)
}

build_proteomics_data <- function( dp, proteomics, gpr )
{
  P = data.frame()
  for(c in unique(dp$condition))
  {
    ### Get last p line for condition c ###
    df = filter(dp, condition==c)
    df = df[dim(df)[1],]
    df = df[,-which(names(df)%in%c("condition", "glc", "phi"))]
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
    #pf$predicted = pf$predicted/sum(pf$predicted)
    pf$predicted = pf$predicted*0.4
    pf$condition = rep(c, dim(pf)[1])
    P            = rbind(P, pf)
  }
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
conditions     = seq(1,30)
GLC            = c(5.0, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0.078125, 0.0390625, 0.01953125, 0.009765625, 0.0048828125, 0.00244140625, 0.001220703125, 0.0006103515625, 0.00030517578125, 0.000152587890625, 7.62939453125e-05, 3.814697265625e-05, 1.9073486328125e-05, 9.5367431640625e-06, 4.76837158203125e-06, 2.384185791015625e-06, 1.1920928955078125e-06, 5.960464477539062e-07, 2.980232238769531e-07, 1.4901161193847656e-07, 7.450580596923828e-08, 3.725290298461914e-08, 1.862645149230957e-08, 9.313225746154785e-09)
mass_fractions = load_mass_fractions()
proteomics     = load_proteomics()
gpr            = load_GPR()

#-------------------------------#
# 2) Load optimum results       #
#-------------------------------#

dstate     = read.table("output/MMSYN_state_optimum.csv", h=T, sep=";")
df         = read.table("output/MMSYN_f_optimum.csv", h=T, sep=";")
dc         = read.table("output/MMSYN_c_optimum.csv", h=T, sep=";")
dv         = read.table("output/MMSYN_v_optimum.csv", h=T, sep=";")
dp         = read.table("output/MMSYN_p_optimum.csv", h=T, sep=";")
db         = read.table("output/MMSYN_b_optimum.csv", h=T, sep=";")
dstate$glc = GLC
df$glc     = GLC
dc$glc     = GLC
dv$glc     = GLC
dp$glc     = GLC
db$glc     = GLC
dp$phi     = dp$Ribosome/rowSums(dp[,-which(names(dp)%in%c("condition", "glc"))])
dp$mu      = dstate$mu
dmf        = build_mass_fraction_data(db, mass_fractions)
dpr        = build_proteomics_data(dp, proteomics, gpr)

p1 = ggplot(dstate, aes(glc, mu)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point() +
  xlab("Glucose (g/L)") +
  ylab("Growth rate") +
  ggtitle("Growth rate depending on glucose concentration") +
  theme_classic()

p2 = ggplot(dp, aes(mu, phi)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, lty=2) +
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2) +
  xlab("Growth rate") +
  ylab("Phi") +
  ggtitle("Growth law") +
  theme_classic()
#plot(db$condition, db$h2o)

p3 = ggplot(filter(dmf, condition==1), aes(observed, predicted)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method="lm", se=F) +
  geom_line(aes(observed, observed), lty=3) +
  geom_line(aes(observed, observedp3), lty=2) +
  geom_line(aes(observed, observedm3), lty=2) +
  geom_point() +
  xlab("Observed mass fractions") +
  ylab("Predicted mass fractions") +
  ggtitle("Observed vs. predicted mass fractions") +
  theme_classic()

p4 = ggplot(filter(dpr, condition==1), aes(observed, predicted)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method="lm", se=F) +
  geom_line(aes(observed, observed), lty=2) +
  geom_line(aes(observed, observedp3), lty=3) +
  geom_line(aes(observed, observedm3), lty=3) +
  geom_point() +
  xlab("Observed proteome fractions") +
  ylab("Predicted proteome fractions") +
  ggtitle("Observed vs. predicted proteome fractions") +

  theme_classic()

plot_grid(p1, p2, p3, p4, ncol=2, labels="AUTO")

