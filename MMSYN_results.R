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
library("ggpmisc")

### Load mass fractions data ###
load_mass_fractions <- function()
{
  d = read.table("/Users/charlesrocabert/git/charlesrocabert/GBA_MMSYN/data/source/Breuer-et-al-2019/mass_fractions.csv", sep=";", h=T, check.names=F)
  rownames(d) = d$ID
  return(d)
}

### Load proteomics data ###
load_proteomics <- function()
{
  d = read.table("/Users/charlesrocabert/git/charlesrocabert/GBA_MMSYN/data/fba_model/MMSYN_proteomics.csv", sep=";", h=T, check.names=F)
  return(d)
}

### Load the GPR ###
load_GPR <- function()
{
  d = read.table("/Users/charlesrocabert/git/charlesrocabert/GBA_MMSYN/data/fba_model/MMSYN_GPR.csv", sep=";", h=T, check.names=F)
  return(d)
}

### Build the mass fraction data ###
build_mass_fraction_data <- function( dc, mass_fractions )
{
  M = data.frame()
  for(c in unique(dc$condition))
  {
    dlocal           = filter(dc, condition==c)
    dlocal           = dlocal[,-which(names(dlocal)%in%c("condition", "glc", "h2o"))]
    dsum             = rowSums(dlocal)
    dlocal           = dlocal/dsum
    dlocal           = dlocal[,which(names(dlocal)%in%mass_fractions$ID)]
    dlocal           = data.frame(names(dlocal), t(dlocal[1,]))
    names(dlocal)    = c("ID", "predicted")
    dlocal$observed  = mass_fractions[dlocal$ID,"Fraction"]/100
    dlocal$condition = rep(c, dim(dlocal)[1])
    dlocal$observed  = dlocal$observed/sum(dlocal$observed)
    dlocal$predicted = dlocal$predicted/sum(dlocal$predicted)
    dlocal$rank      = rank(dlocal$observed)
    M                = rbind(M, dlocal)
  }
  M$condition   = factor(M$condition, levels=unique(sort(as.integer(M$condition))))
  M$observedp3  = M$observed*3
  M$observedm3  = M$observed/3
  M$observedp10 = M$observed*10
  M$observedm10 = M$observed/10
  return(M)
}

build_proteomics_data <- function( dp, proteomics, gpr, proteome_fractions, model_name )
{
  P = data.frame()
  for(c in unique(dp$condition))
  {
    ### Get last p line for condition c ###
    dlocal = filter(dp, condition==c)
    dlocal = dlocal[,-which(names(dlocal)%in%c("condition", "glc", "phi", "mu"))]
    #df  = df#/sum(df)
    ### Parse the GPR to build the proteome masses ###
    plocal           = data.frame(proteomics$protein, rep(0.0, length(proteomics$protein)))
    names(plocal)    = c("protein", "predicted")
    plocal$observed  = proteomics$mass*0.001#/sum(proteomics$mass)
    rownames(plocal) = plocal$protein
    for(i in seq(1, dim(gpr)[1]))
    {
      reaction    = gpr$reaction[i]
      step1       = gpr$GPR[i]
      step2       = strsplit(step1,"|",fixed=T)
      total_stoic = 0.0
      #print(paste("> Reaction ", reaction, " --> gpr =", step1))
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
        if (reaction%in%names(dlocal))
        {
          plocal[protein,"predicted"] = plocal[protein,"predicted"]+dlocal[,reaction]*stoic/total_stoic
        }
      }
    }
    plocal          = filter(plocal, predicted > 0.0)
    plocal          = filter(plocal, observed > 0.0)
    #plocal$observed  = plocal$observed/sum(plocal$observed)
    plocal$predicted = plocal$predicted/sum(plocal$predicted)
    plocal$predicted = plocal$predicted*proteome_fractions[[model_name]]
    plocal$condition = rep(c, dim(plocal)[1])
    P                = rbind(P, plocal)
  }
  P$condition   = factor(P$condition, levels=unique(sort(as.integer(P$condition))))
  P$observedp3  = P$observed*3
  P$observedm3  = P$observed/3
  P$observedp10 = P$observed*10
  P$observedm10 = P$observed/10
  return(P)
}

plot_trends_for_given_model <- function( model_name )
{
  dstate = read.table(paste0("output/",model_name,"_state_optimum.csv"), h=T, sep=";")
  df     = read.table(paste0("output/",model_name,"_f_optimum.csv"), h=T, sep=";")
  db     = read.table(paste0("output/",model_name,"_b_optimum.csv"), h=T, sep=";")
  dc     = read.table(paste0("output/",model_name,"_c_optimum.csv"), h=T, sep=";")
  dv     = read.table(paste0("output/",model_name,"_v_optimum.csv"), h=T, sep=";")
  dp     = read.table(paste0("output/",model_name,"_p_optimum.csv"), h=T, sep=";")
  ##########################################
  # dstate$glc = GLC
  # df$glc     = GLC
  # dc$glc     = GLC
  # dv$glc     = GLC
  # dp$glc     = GLC
  # db$glc     = GLC
  dp$phi     = dp$Ribosome/rowSums(dp[,-which(names(dp)%in%c("condition", "glc"))])
  dp$mu      = dstate$mu
  dmf        = build_mass_fraction_data(dc, mass_fractions)
  dpr        = build_proteomics_data(dp, proteomics, gpr, proteome_fractions, model_name)
  p = ggplot(dmf, aes(observed, predicted, color=condition)) +
    scale_x_log10() + scale_y_log10() +
    geom_abline(intercept=0, slope=1) +
    geom_point() +
    geom_smooth(method="lm", se=F) +
    facet_wrap(.~condition)
  return(p)
  ##########################################
  p1 = ggplot(dstate, aes(condition, mu)) +
    geom_hline(aes(yintercept=0.4, color="Experimental"), lty=2) +
    geom_hline(aes(yintercept=0.34, color="FBA"), lty=2) +
    #scale_x_log10() +
    #scale_y_log10() +
    geom_point() +
    xlab("Glucose (g/L)") +
    ylab("Growth rate (h-1)") +
    ggtitle("Growth rate VS. condition") +
    theme_classic() +
    theme(legend.position="bottom")
  ##########################################
  p2 = ggplot(dp, aes(mu, phi)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, lty=2) +
    #geom_hline(yintercept=0, lty=2) +
    #geom_vline(xintercept=0, lty=2) +
    xlab("Growth rate (h-1)") +
    ylab("Phi (ribosome proteome fraction)") +
    ggtitle("Growth law") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    theme_classic()
  ##########################################
  p3 = ggplot(filter(dmf, condition==1), aes(observed, predicted)) +
    scale_x_log10() +
    scale_y_log10() +
    #geom_smooth(method="lm", se=F) +
    stat_poly_line() +
    stat_poly_eq(use_label("R2", "P")) +
    geom_line(aes(observed, observed), lty=3) +
    geom_line(aes(observed, observedp3), lty=2) +
    geom_line(aes(observed, observedm3), lty=2) +
    geom_point() +
    xlab("Observed mass fractions (% of gDW)") +
    ylab("Predicted mass fractions (% of gDW)") +
    ggtitle("Observed vs. predicted mass fractions") +
    theme_classic()
  ##########################################
  p4 = ggplot(filter(dpr, condition==1), aes(observed, predicted)) +
    scale_x_log10() +
    scale_y_log10() +
    #geom_smooth(method="lm", se=F) +
    stat_poly_line() +
    stat_poly_eq(use_label("R2", "P")) +
    geom_line(aes(observed, observed), lty=2) +
    geom_line(aes(observed, observedp3), lty=3) +
    geom_line(aes(observed, observedm3), lty=3) +
    geom_point() +
    xlab("Observed proteome fractions") +
    ylab("Predicted proteome fractions\n(weighted by dumped fraction)") +
    ggtitle("Observed vs. predicted proteome fractions") +
    theme_classic()
  ##########################################
  p = plot_grid(p1, p2, p3, p4, ncol=2, labels="AUTO")
  return(p)
}

plot_growth_law <- function( model_name )
{
  dstate = read.table(paste0("output/",model_name,"_state_optimum.csv"), h=T, sep=";")
  dp     = read.table(paste0("output/",model_name,"_p_optimum.csv"), h=T, sep=";")
  ##########################################
  dp$glc = GLC
  dp$phi = dp$Ribosome/rowSums(dp[,-which(names(dp)%in%c("condition", "glc"))])
  dp$mu  = dstate$mu
  ##########################################
  p = ggplot(dp, aes(mu, phi)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, lty=2) +
    #geom_hline(yintercept=0, lty=2) +
    #geom_vline(xintercept=0, lty=2) +
    xlab("Growth rate (h-1)") +
    ylab("Phi (ribosome proteome fraction)") +
    ggtitle(paste0("Growth law (model ",model_name,")")) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    theme_classic()
  ##########################################
  return(p)
}

compare_growth_rates <- function( models, selected_condition )
{
  D = data.frame()
  for (model_name in models)
  {
    dstate = read.table(paste0("output/",model_name,"_state_optimum.csv"), h=T, sep=";")
    D = rbind(D, filter(dstate, condition==selected_condition))
  }
  D$model = models
  return(D)
}

compare_mass_fraction_predictions <- function( models, selected_condition )
{
  D = data.frame()
  for (model_name in models)
  {
    dc  = read.table(paste0("output/",model_name,"_c_optimum.csv"), h=T, sep=";")
    dmf = build_mass_fraction_data(dc, mass_fractions)
    #print(head(dmf))
    reg  = lm(log10(dmf$predicted)~log10(dmf$observed))
    r2   = summary(reg)$adj.r.squared
    pval = summary(reg)$coefficients[2,4]
    #print(pval)
    D = rbind(D, c(r2, pval))
  }
  names(D) = c("r2", "pvalue")
  D$model = models
  return(D)
}

compare_proteome_fraction_predictions <- function( models, selected_condition )
{
  D = data.frame()
  for (model_name in models)
  {
    dp  = read.table(paste0("output/",model_name,"_p_optimum.csv"), h=T, sep=";")
    dpr = build_proteomics_data(dp, proteomics, gpr, proteome_fractions, model_name)
    #print(head(dmf))
    reg  = lm(log10(dpr$predicted)~log10(dpr$observed))
    r2   = summary(reg)$adj.r.squared
    pval = summary(reg)$coefficients[2,4]
    #print(pval)
    D = rbind(D, c(r2, pval))
  }
  names(D) = c("r2", "pvalue")
  D$model = models
  return(D)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/charlesrocabert/GBA_Evolution_2/")

#-------------------------------#
# 1) Load experimental datasets #
#-------------------------------#
proteome_fractions = list(
  "MMSYN_1111" = 0.33842338423384233,
  "MMSYN_1110" = 0.33842338423384233,
  "MMSYN_1101" = 0.33842338423384233,
  "MMSYN_1100" = 0.33842338423384233,
  "MMSYN_1011" = 1.0,
  "MMSYN_1010" = 1.0,
  "MMSYN_1001" = 1.0,
  "MMSYN_1000" = 1.0,
  "MMSYN_0111" = NA,
  "MMSYN_0110" = 0.35245352453524537,
  "MMSYN_0101" = NA,
  "MMSYN_0100" = 0.35245352453524537,
  "MMSYN_0011" = NA,
  "MMSYN_0010" = 1.0,
  "MMSYN_0001" = NA,
  "MMSYN_0000" = 1.0
)
conditions     = seq(1,30)
GLC            = c(5.0, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0.078125, 0.0390625, 0.01953125, 0.009765625, 0.0048828125, 0.00244140625, 0.001220703125, 0.0006103515625, 0.00030517578125, 0.000152587890625, 7.62939453125e-05, 3.814697265625e-05, 1.9073486328125e-05, 9.5367431640625e-06, 4.76837158203125e-06, 2.384185791015625e-06, 1.1920928955078125e-06, 5.960464477539062e-07, 2.980232238769531e-07, 1.4901161193847656e-07, 7.450580596923828e-08, 3.725290298461914e-08, 1.862645149230957e-08, 9.313225746154785e-09)
mass_fractions = load_mass_fractions()
proteomics     = load_proteomics()
gpr            = load_GPR()
model_name     = "MMSYN_1010"
models         = c("MMSYN_0000", "MMSYN_1000", "MMSYN_0100", "MMSYN_0010",
                   "MMSYN_0110", "MMSYN_1100", "MMSYN_1010", "MMSYN_1110")


#-------------------------------#
# 2) Load results per model     #
#-------------------------------#
# plot_trends_for_given_model("MMSYN_0000")
# plot_trends_for_given_model("MMSYN_1000")
# plot_trends_for_given_model("MMSYN_0100")
# plot_trends_for_given_model("MMSYN_0010")
# plot_trends_for_given_model("MMSYN_0110")
# plot_trends_for_given_model("MMSYN_1100")
# plot_trends_for_given_model("MMSYN_1010")
# plot_trends_for_given_model("MMSYN_1110")

plot_trends_for_given_model("MMSYN_1011")

stop()

p1 = plot_growth_law("MMSYN_0000")
p2 = plot_growth_law("MMSYN_1000")
p3 = plot_growth_law("MMSYN_0100")
p4 = plot_growth_law("MMSYN_0010")
p5 = plot_growth_law("MMSYN_0110")
p6 = plot_growth_law("MMSYN_1100")
p7 = plot_growth_law("MMSYN_1010")
p8 = plot_growth_law("MMSYN_1110")
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, labels="AUTO", ncol=3)

########################################################
########################################################
########################################################

D1 = compare_growth_rates(models, "1")
D1 = D1[order(D1$mu, decreasing=T),]
D1$model = factor(D1$model, levels=D1$model)
ggplot(D1 ,aes(model, mu)) +
  geom_hline(aes(yintercept=0.4, color="Experimental"), lty=2) +
  geom_hline(aes(yintercept=0.34, color="FBA"), lty=2) +
  geom_point() +
  xlab("") +
  ylab("Growth rate (h-1)") +
  ggtitle("Growth rate with 5g/L glucose") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))

D2 = compare_mass_fraction_predictions(models, "1")
D2 = D2[order(D2$r2, decreasing=T),]
D2$model = factor(D2$model, levels=D2$model)
p1 = ggplot(D2 ,aes(model, r2)) +
  geom_point() +
  xlab("") +
  ylab("Linear regression R^2") +
  ggtitle("Predicted VS. observed mass fractions with 5g/L glucose") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
p2 = ggplot(D2 ,aes(model, log10(pvalue))) +
  geom_point() +
  xlab("") +
  ylab("Associated p-value") +
  ggtitle("Predicted VS. observed mass fractions with 5g/L glucose") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
plot_grid(p1, p2, labels="AUTO", ncol=2)

D3 = compare_proteome_fraction_predictions(models, "1")
D3 = D3[order(D3$r2, decreasing=T),]
D3$model = factor(D3$model, levels=D3$model)
p1 = ggplot(D3 ,aes(model, r2)) +
  geom_point() +
  xlab("") +
  ylab("Linear regression R^2") +
  ggtitle("Predicted VS. observed proteome fractions\nwith 5g/L glucose") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
p2 = ggplot(D3 ,aes(model, log10(pvalue))) +
  geom_point() +
  xlab("") +
  ylab("Associated p-value") +
  ggtitle("Predicted VS. observed proteome fractions\nwith 5g/L glucose") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
plot_grid(p1, p2, labels="AUTO", ncol=2)

