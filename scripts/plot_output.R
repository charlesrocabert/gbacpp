#!/usr/bin/Rscript
# coding: utf-8

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

### Load the experimental dataset ###
d_obs = read.table("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/wet_experiments/observed_mu.csv", sep=";", h=T)
Davg = tapply(d_obs$mu, d_obs$glc, mean)
Davg = data.frame(glc=as.numeric(names(Davg)), mu=Davg)

### Build the glucose vector ###
glc_obs = c(0.0, 1e-05, 0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0)
glc_obs = glc_obs+0.4
glc_sim = 10.0^seq(-5.9, 1, by=0.1)
glc_vec = c(glc_obs, glc_sim)


d_state = read.table("../output/mmsyn_fcr_state_optimum.csv", h=T, sep=";")
d_b     = read.table("../output/mmsyn_fcr_b_optimum.csv", h=T, sep=";")

d_state$glc = glc_vec
d_b$glc     = glc_vec

MF_cor = data.frame()
for (i in seq(1, dim(d_b)[1]))
{
  res     = mass_fractions_cor(d_b, i )
  MF_cor = rbind(MF_cor, res)
}
MF_cor$condition = d_b$condition
MF_cor$glc       = glc_vec
names(MF_cor) = c("r2", "pval", "R2", "condition", "glc")

ggplot(d_state, aes(glc, mu)) +
  geom_hline(yintercept=mean(d_obs$mu), col="pink") +
  geom_line() +
  geom_point(data=Davg, aes(glc, mu, color="Observed")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic()

mean(d_obs$mu)
median(d_obs$mu)
max(d_state$mu)
