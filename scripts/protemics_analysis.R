#!/usr/bin/Rscript
# coding: utf-8

library("tidyverse")
library("rstudioapi")
library("cowplot")
library("ggpmisc")
library("Matrix")
library("ggrepel")
library("latex2exp")

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
  names(res)    = c("p_id", "mass")
  rownames(res) = res$p_id
  for (r_id in r_ids)
  {
    e_conc  = X[r_id][[1]]
    contrib = protein_contributions[protein_contributions$reaction==r_id,]
    if (dim(contrib)[1] > 0)
    {
      for (j in seq(1, dim(contrib)[1]))
      {
        p_id             = contrib$protein[j]
        mass             = e_conc*contrib$contribution[j]
        res[p_id,"mass"] = res[p_id,"mass"]+mass
      }
    }
  }
  res$fraction = res$mass/sum(res$mass)
  return(res)
}

### Collect proteins excluded from the model ###
collect_excluded_proteins <- function( list_of_reactions, enzyme_composition )
{
  excluded_proteins = c()
  for(i in seq(1, dim(enzyme_composition)[1]))
  {
    r_id = enzyme_composition$info_sample_rid[i]
    if (!r_id %in% list_of_reactions)# | r_id %in% c("Ribosome", "AAabc", "AATRS"))
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
  Xsum = sum(obs_proteomics$mass)
  Xsim = sum(filter(obs_proteomics, protein%in%rownames(sim_proteomics))$mass)
  print(paste("> Modeled proteome fraction", Xsim/Xsum))
  D              = data.frame(sim_proteomics$p_id, sim_proteomics$mass, sim_proteomics$fraction)
  names(D)       = c("id", "sim_mass", "sim_fraction")
  D              = filter(D, id%in%obs_proteomics$protein)
  D$obs_mass     = obs_proteomics[D$id, "mass"]
  D              = D[!D$id%in%excluded,]
  D$obs_fraction = D$obs_mass/sum(D$obs_mass)
  D$sim_fraction = D$sim_mass/sum(D$sim_mass)#*0.3435353
  D$obsp3        = D$obs_fraction*3
  D$obsm3        = D$obs_fraction/3
  D$obsp10       = D$obs_fraction*10
  D$obsm10       = D$obs_fraction/10
  return(D)
}

### Correlation of proteomics prediction ###
proteomics_cor <- function( d_p )
{
  D     = build_proteomics_data(d_p, dim(d_p)[1])
  reg   = lm(log10(D$sim_fraction)~log10(D$obs_fraction))
  r2    = summary(reg)$adj.r.squared
  pval  = summary(reg)$coefficients[,4][[2]]
  SSres = sum((log10(D$sim_fraction)-log10(D$obs_fraction))^2)
  M     = mean(log10(D$obs_fraction))
  SStot = sum((log10(D$obs_fraction)-M)^2)
  R2    = 1-SSres/SStot
  return(c(r2, pval, R2))
}

### Evolution of proteomics prediction ###
proteomics_evolution <- function( d_p, step )
{
  iter     = c()
  cor_vec  = c()
  pval_vec = c()
  R2_vec   = c()
  for(i in seq(1, dim(d_p)[1], by=step))
  {
    iter     = c(iter, i)
    D        = build_proteomics_data(d_p, i)
    res      = proteomics_cor(D)
    cor_vec  = c(cor_vec, res[1])
    pval_vec = c(pval_vec, res[2])
    R2_vec   = c(R2_vec, res[3])
  }
  D = data.frame(iter, cor_vec, pval_vec, R2_vec)
  names(D) = c("index", "r2", "pval", "R2")
  return(D)
}

### Plot predicted proteomics ###
plot_predicted_proteomics <- function( pr_data, R2, r2 )
{
  p = ggplot(pr_data, aes(obs_fraction, sim_fraction)) +
    # geom_abline(slope=1, intercept=0, color="pink") +
    # geom_line(aes(obs_fraction, obsp3), color="grey", lty=2) +
    # geom_line(aes(obs_fraction, obsm3), color="grey", lty=2) +
    # geom_line(aes(obs_fraction, obsp10), color="grey", lty=3) +
    # geom_line(aes(obs_fraction, obsm10), color="grey", lty=3) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_x_log10() + scale_y_log10() +
    annotate("text", x=0.001, y=0.05, label=paste0("italic(R)^2", "==", round(R2,5)), hjust=0, parse=T) +
    annotate("text", x=0.001, y=0.04, label=paste0("italic(r)^2", "==", round(r2,5)), hjust=0, parse=T) +
    #geom_text_repel(aes(label=id), size = 3.5) +
    xlab("Observed") +
    ylab("Simulated") +
    ggtitle("Proteomics") +
    theme_classic()
  return(p)
}

### Plot proteomics correlation during optimization ###
plot_proteomics_evolution <- function( pr_evol_data )
{
  p1 = ggplot(pr_evol_data, aes(index, r2)) +
    geom_line() +
    xlab("Iteration") +
    ylab("Linear regression adj. R-squared") +
    ggtitle("Correlation with observed proteomics") +
    theme_classic()
  p2 = ggplot(pr_evol_data, aes(index, pval)) +
    geom_line() +
    geom_hline(yintercept=0.05, col="pink") +
    scale_y_log10() +
    xlab("Iteration") +
    ylab("Linear regression p-value") +
    ggtitle("Correlation with observed proteomics") +
    theme_classic()
  p3 = ggplot(pr_evol_data, aes(index, R2)) +
    geom_line() +
    xlab("Iteration") +
    ylab("R2") +
    ggtitle("R2 with observed proteomics") +
    theme_classic()
  return(list(p1, p2, p3))
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
mu_obs       = c(0.4185156420108013, 0.34575412067094446, 0.32362916852088447, 0.28943509125099154, 0.30019744370357504, 0.358911043059261, 0.3009631588728481, 0.42854660068743594, 0.29781579360835236, 0.3927863575066446, 0.3775322888686228, 0.4893620794586438, 0.47200997079180596, 0.5141549673685587, 0.5092329943189714, 0.497863115357777, 0.5181987729683097, 0.5010039120856962, 0.46432998072355164, 0.43444304843015225)
cond         = seq(1,20)
d_obs        = data.frame(cond, x_adjust, mu_obs)
names(d_obs) = c("condition", "glc", "mu")

model_path  = "../csv_models"
model_path  = "../../gbapy/tutorials/MMSYN_tutorial/models"
output_path = "../hpc_download"
output_path = "/Users/charlesrocabert/Desktop/mmsyn_output"
model_name  = "mmsyn_fcr"
condition   = 20

d_state = read.table(paste0(output_path,"/",model_name,"_",condition,"_state_trajectory.csv"), h=T, sep=";", check.names=F)
d_f     = read.table(paste0(output_path,"/",model_name,"_",condition,"_f_trajectory.csv"), h=T, sep=";", check.names=F)
d_v     = read.table(paste0(output_path,"/",model_name,"_",condition,"_v_trajectory.csv"), h=T, sep=";", check.names=F)
d_p     = read.table(paste0(output_path,"/",model_name,"_",condition,"_p_trajectory.csv"), h=T, sep=";", check.names=F)
d_b     = read.table(paste0(output_path,"/",model_name,"_",condition,"_b_trajectory.csv"), h=T, sep=";", check.names=F)
d_c     = read.table(paste0(output_path,"/",model_name,"_",condition,"_c_trajectory.csv"), h=T, sep=";", check.names=F)

MF      = build_mass_fractions_data(d_b, dim(d_b)[1])
MF_cor  = mass_fractions_cor(d_b)
MF_evol = mass_fraction_evolution(d_b, 10)
PR      = build_proteomics_data(d_p, dim(d_p)[1])
PR_cor  = proteomics_cor(d_p)

p1 = plot_growth_rate(d_state, d_obs)
p2 = plot_protein_fraction(d_c)
p3 = plot_mass_fractions(MF, MF_cor[3], MF_cor[1])
p4 = plot_proteomics(PR, PR_cor[3], PR_cor[1])
p_mf = plot_mass_fractions_evolution(MF_evol)
plot_grid(p1, p2, p3, p4, p_mf[[1]], p_mf[[3]], ncol=2)

X = d_f[dim(d_f)[1],-which(names(d_f)%in%c("condition", "iter", "t", "dt"))]
X = data.frame(names(X), t(X))
names(X) = c("reaction_id", "flux_fraction")
X = X[order(X[,2], decreasing=T),]
tail(X)

