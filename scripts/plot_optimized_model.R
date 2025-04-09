#!/usr/bin/Rscript
# coding: utf-8

library("tidyverse")
library("rstudioapi")
library("cowplot")
library("ggpmisc")
library("Matrix")
library("ggrepel")
library("latex2exp")
library("MASS")

### Load the observed mass fractions data ###
load_observed_mass_fractions <- function()
{
  d           = read.table("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/manual_curation/MMSYN_mass_fractions.csv", sep=";", h=T, check.names=F)
  rownames(d) = d$ID
  return(d)
}

### Build the mass fraction data ###
build_mass_fractions_data <- function( d_b, i )
{
  d_mf           = load_observed_mass_fractions()
  D              = d_b[i,-which(names(d_b)%in%c("condition", "iter", "t", "dt", "h2o"))]
  D              = D[,which(names(D)%in%d_mf$ID)]
  D              = data.frame(names(D), t(D))
  names(D)       = c("id", "sim_mass")
  D$obs_fraction = d_mf[D$id,"Fraction"]
  D$obs_fraction = D$obs*0.01
  D$sim_fraction = D$sim_mass/sum(D$sim_mass)
  D$obsp3        = D$obs_fraction*3
  D$obsm3        = D$obs_fraction/3
  D$obsp10       = D$obs_fraction*10
  D$obsm10       = D$obs_fraction/10
  return(D)
}

### Correlation of mass fraction prediction ###
mass_fractions_correlation <- function( d_b, i )
{
  D     = build_mass_fractions_data(d_b, i)
  reg   = lm(log10(D$sim_fraction)~log10(D$obs_fraction))
  r2    = summary(reg)$adj.r.squared
  pval  = summary(reg)$coefficients[,4][[2]]
  SSres = sum((log10(D$sim_fraction)-log10(D$obs_fraction))^2)
  M     = mean(log10(D$obs_fraction))
  SStot = sum((log10(D$obs_fraction)-M)^2)
  R2    = 1-SSres/SStot
  return(c(r2, pval, R2))
}

### Evolution of mass fraction prediction ###
mass_fractions_optimization <- function( d_b, step )
{
  iter     = c()
  cor_vec  = c()
  pval_vec = c()
  R2_vec   = c()
  for(i in seq(1, dim(d_b)[1], by=step))
  {
    iter     = c(iter, i)
    res      = mass_fractions_correlation(d_b, i)
    cor_vec  = c(cor_vec, res[1])
    pval_vec = c(pval_vec, res[2])
    R2_vec   = c(R2_vec, res[3])
  }
  if (iter[length(iter)] != dim(d_b)[1])
  {
    iter     = c(iter, dim(d_b)[1])
    res      = mass_fractions_correlation(d_b, dim(d_b)[1])
    cor_vec  = c(cor_vec, res[1])
    pval_vec = c(pval_vec, res[2])
    R2_vec   = c(R2_vec, res[3])
  }
  D = data.frame(iter, cor_vec, pval_vec, R2_vec)
  names(D) = c("index", "r2", "pval", "R2")
  return(D)
}

### Load observed proteomics data ###
load_observed_proteomics <- function()
{
  d           = read.table(paste0("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/manual_curation/MMSYN_proteomics.csv"), sep=";", h=T, check.names=F)
  d$protein   = str_replace(d$locus, "JCVISYN3A", "protein")
  d$obs_mass  = d$mg_per_gP*0.001
  rownames(d) = d$protein
  d           = filter(d, obs_mass > 0.0)
  return(d)
}

### Load protein contributions ###
load_protein_contributions <- function( model_path, model_name )
{
  filename    = paste0(model_path,"/",model_name,"/protein_contributions.csv")
  d           = read.table(filename, h=T, sep=";", check.names=F)
  nrow        = length(unique(d$protein))
  ncol        = length(unique(d$reaction))
  M           = matrix(rep(0.0, nrow*ncol), nrow, ncol)
  rownames(M) = sort(unique(d$protein))
  colnames(M) = sort(unique(d$reaction))
  for(i in seq(1, dim(d)[1]))
  {
    M[d$protein[i], d$reaction[i]] = d$contribution[i]
  }
  return(M)
}

### Load enzyme composition ###
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
calculate_simulated_proteomics <- function( model_path, model_name, d_p, i )
{
  X                = d_p[i, -which(names(d_p)%in%c("condition", "iter", "t", "dt"))]
  M                = load_protein_contributions(model_path, model_name)
  X                = as.vector(t(X))
  Y                = M%*%X
  r_ids            = colnames(M)
  p_ids            = rownames(M)
  res              = data.frame(p_ids, Y)
  names(res)       = c("p_id", "sim_mass")
  rownames(res)    = res$p_id
  res$sim_fraction = res$sim_mass/sum(res$sim_mass)
  return(res)
}

### Collect proteins excluded from the model ###
collect_excluded_proteins <- function( list_of_reactions, enzyme_composition )
{
  excluded_proteins = c()
  for(i in seq(1, dim(enzyme_composition)[1]))
  {
    r_id = enzyme_composition$info_sample_rid[i]
    if (!r_id %in% list_of_reactions)# | r_id %in% c("Ribosome"))#, "AAabc", "AATRS"))
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
build_proteomics_data <- function( model_path, model_name, d_p, i )
{
  obs_proteomics = load_observed_proteomics()
  enz_comp       = load_enzyme_composition()
  reaction_ids   = names(d_p[,-which(names(d_p)%in%c("condition", "iter", "t", "dt"))])
  excluded       = collect_excluded_proteins(reaction_ids, enz_comp)
  sim_proteomics = calculate_simulated_proteomics(model_path, model_name, d_p, i)
  Xsum           = sum(obs_proteomics$obs_mass)
  Xsim           = sum(filter(obs_proteomics, protein%in%rownames(sim_proteomics))$obs_mass)
  print(paste("> Modeled proteome fraction", Xsim/Xsum))
  D              = data.frame(sim_proteomics$p_id, sim_proteomics$sim_mass, sim_proteomics$sim_fraction)
  names(D)       = c("id", "sim_mass", "sim_fraction")
  D              = filter(D, id%in%obs_proteomics$protein)
  D$obs_mass     = obs_proteomics[D$id, "obs_mass"]
  D              = D[!D$id%in%excluded,]
  D$obs_fraction = D$obs_mass/sum(D$obs_mass)
  D$sim_fraction = D$sim_mass/sum(D$sim_mass)
  D$obsp3        = D$obs_fraction*3
  D$obsm3        = D$obs_fraction/3
  D$obsp10       = D$obs_fraction*10
  D$obsm10       = D$obs_fraction/10
  return(D)
}

### Correlation of proteomics prediction ###
proteomics_correlation <- function( model_path, model_name, d_p, i )
{
  D     = build_proteomics_data(model_path, model_name, d_p, i)
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
proteomics_optimization <- function( model_path, model_name, d_p, step )
{
  iter     = c()
  cor_vec  = c()
  pval_vec = c()
  R2_vec   = c()
  for(i in seq(1, dim(d_p)[1], by=step))
  {
    iter     = c(iter, i)
    res      = proteomics_correlation(model_path, model_name, d_p, i)
    cor_vec  = c(cor_vec, res[1])
    pval_vec = c(pval_vec, res[2])
    R2_vec   = c(R2_vec, res[3])
  }
  if (iter[length(iter)] != dim(d_p)[1])
  {
    iter     = c(iter, dim(d_p)[1])
    res      = proteomics_correlation(model_path, model_name, d_p, dim(d_p)[1])
    cor_vec  = c(cor_vec, res[1])
    pval_vec = c(pval_vec, res[2])
    R2_vec   = c(R2_vec, res[3])
  }
  D = data.frame(iter, cor_vec, pval_vec, R2_vec)
  names(D) = c("index", "r2", "pval", "R2")
  return(D)
}

### Plot predicted growth rates ###
plot_predicted_growth_rate <- function( d_state, d_obs )
{
  last_mu = d_state$mu[dim(d_state)[1]]
  p = ggplot(d_state, aes(iter, mu)) +
    geom_line() +
    geom_hline(yintercept=max(d_obs$mu), color="pink") +
    xlab("Iterations") +
    ylab("Growth rate") +
    ggtitle(paste0("Growth rate (\u03BC = ",round(last_mu,3),")")) +
    theme_classic()
  return(p)
}

### Plot predicted protein fraction ###
plot_predicted_protein_fraction <- function( d_c )
{
  dl                  = d_c[,-which(names(d_c)%in%c("condition","iter","t","dt", "h2o"))]
  dl$sum              = rowSums(dl)
  dl$Protein_fraction = dl$Protein/dl$sum
  dl$iter             = d_c$iter
  last_prot_fraction  = dl$Protein_fraction[dim(dl)[1]]
  p = ggplot(dl, aes(iter, Protein_fraction)) +
    geom_hline(yintercept=0.54727, color="pink") +
    geom_line() +
    xlab("Iterations") +
    ylab("Protein fraction") +
    ggtitle(paste0("Protein fraction (Pf = ",round(last_prot_fraction,3),", obs = ",0.547,")")) +
    ylim(0.5,1) +
    theme_classic()
  return(p)
}

### Plot predicted mass fractions ###
plot_predicted_mass_fractions <- function( mf_data, R2, r2 )
{
  p = ggplot(mf_data, aes(obs_fraction, sim_fraction)) +
    geom_abline(slope=1, intercept=0, color="pink") +
    geom_line(aes(obs_fraction, obsp3), color="grey", lty=2) +
    geom_line(aes(obs_fraction, obsm3), color="grey", lty=2) +
    geom_line(aes(obs_fraction, obsp10), color="grey", lty=3) +
    geom_line(aes(obs_fraction, obsm10), color="grey", lty=3) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_x_log10() + scale_y_log10() +
    #geom_text_repel(aes(label=id), size = 3.5) +
    annotate("text", x=1e-5, y=5e-1, label=paste0("italic(R)^2", "==", round(R2,5)), hjust=0, parse=T) +
    annotate("text", x=1e-5, y=5e-2, label=paste0("italic(r)^2", "==", round(r2,5)), hjust=0, parse=T) +
    xlab("Observed") +
    ylab("Simulated") +
    ggtitle("Metabolite mass fractions") +
    theme_classic()
  return(p)
}

### Plot mass fraction correlation during optimization ###
plot_mass_fractions_optimization <- function( mf_evol_data )
{
  p1 = ggplot(mf_evol_data, aes(index, r2)) +
    geom_line() +
    xlab("Iteration") +
    ylab("Linear regression adj. R-squared") +
    ggtitle("Correlation with observed mass fractions") +
    theme_classic()
  p2 = ggplot(mf_evol_data, aes(index, pval)) +
    geom_line() +
    geom_hline(yintercept=0.05, col="pink") +
    scale_y_log10() +
    xlab("Iteration") +
    ylab("Linear regression p-value") +
    ggtitle("Correlation with observed mass fractions") +
    theme_classic()
  p3 = ggplot(mf_evol_data, aes(index, R2)) +
    geom_line() +
    xlab("Iteration") +
    ylab("R2") +
    ggtitle("R2 with observed mass fractions") +
    theme_classic()
  return(list(p1, p2, p3))
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
plot_proteomics_optimization <- function( pr_evol_data )
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

### Build experimental conditions ###
build_experimental_conditions <- function()
{
  mu_obs       = c(0.4185156420108013, 0.34575412067094446, 0.32362916852088447, 0.28943509125099154, 0.30019744370357504, 0.358911043059261, 0.3009631588728481, 0.42854660068743594, 0.29781579360835236, 0.3927863575066446, 0.3775322888686228, 0.4893620794586438, 0.47200997079180596, 0.5141549673685587, 0.5092329943189714, 0.497863115357777, 0.5181987729683097, 0.5010039120856962, 0.46432998072355164, 0.43444304843015225)
  glc_obs      = c(0.0, 1e-05, 0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0)
  glc_obs      = glc_obs+0.4
  d_obs        = data.frame(glc_obs+0.4, mu_obs)
  names(d_obs) = c("glc", "mu")
  return(d_obs)
}


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)

### Build the experimental dataset ###
d_obs   = build_experimental_conditions()
glc_sim = 10.0^seq(-5.9, 1, by=0.1)
glc_vec = c(d_obs$glc, glc_sim)

model_path  = "../csv_models"
model_path  = "../../gbapy/tutorials/MMSYN_tutorial/models"
output_path = "../hpc_download"
#output_path = "/Users/charlesrocabert/Desktop/mmsyn_output/v1"
output_path = "/Users/charlesrocabert/Desktop/mmsyn_output/v2_no_transporter_kcat"
model_name  = "mmsyn_fcr"
condition   = 20

d_state = read.table(paste0(output_path,"/",model_name,"_",condition,"_state_trajectory.csv"), h=T, sep=";", check.names=F)
d_f     = read.table(paste0(output_path,"/",model_name,"_",condition,"_f_trajectory.csv"), h=T, sep=";", check.names=F)
d_v     = read.table(paste0(output_path,"/",model_name,"_",condition,"_v_trajectory.csv"), h=T, sep=";", check.names=F)
d_p     = read.table(paste0(output_path,"/",model_name,"_",condition,"_p_trajectory.csv"), h=T, sep=";", check.names=F)
d_b     = read.table(paste0(output_path,"/",model_name,"_",condition,"_b_trajectory.csv"), h=T, sep=";", check.names=F)
d_c     = read.table(paste0(output_path,"/",model_name,"_",condition,"_c_trajectory.csv"), h=T, sep=";", check.names=F)

MF       = build_mass_fractions_data(d_b, dim(d_b)[1])
MF_optim = mass_fractions_optimization(d_b, 100)
PR       = build_proteomics_data(model_path, model_name, d_p, dim(d_p)[1])
PR_optim = proteomics_optimization(model_path, model_name, d_p, 100)

p1 = plot_predicted_growth_rate(d_state, d_obs)
p2 = plot_predicted_protein_fraction(d_c)
p3 = plot_predicted_mass_fractions(MF, MF_optim$R2[dim(MF_optim)[1]], MF_optim$r2[dim(MF_optim)[1]])
p4 = plot_predicted_proteomics(PR, PR_optim$R2[dim(PR_optim)[1]], PR_optim$r2[dim(PR_optim)[1]])
p_mf = plot_mass_fractions_optimization(MF_optim)
plot_grid(p1, p2, p3, p4, p_mf[[1]], p_mf[[3]], ncol=2)

X        = d_f[dim(d_f)[1],-which(names(d_f)%in%c("condition", "iter", "t", "dt"))]
X        = data.frame(names(X), t(X))
names(X) = c("reaction_id", "flux_fraction")
X        = X[order(X[,2], decreasing=T),]

### PHI (Ribosome proteome fraction) ###
Y     = d_p[,-which(names(d_p)%in%c("condition", "iter", "t", "dt"))]
Ysum  = rowSums(Y)
Y$Phi = Y$Ribosome/Ysum
#plot(Y$Phi, type="l")

rib_prot = c("protein_0025", "protein_0027", "protein_0082", "protein_0137", "protein_0148", "protein_0149", "protein_0198", "protein_0199", "protein_0238", "protein_0294", "protein_0362", "protein_0365", "protein_0422", "protein_0482", "protein_0499", "protein_0500", "protein_0501", "protein_0930", "protein_0526", "protein_0540", "protein_0637", "protein_0638", "protein_0644", "protein_0646", "protein_0647", "protein_0648", "protein_0653", "protein_0654", "protein_0655", "protein_0656", "protein_0657", "protein_0658", "protein_0659", "protein_0660", "protein_0661", "protein_0662", "protein_0663", "protein_0664", "protein_0665", "protein_0666", "protein_0667", "protein_0668", "protein_0669", "protein_0670", "protein_0671", "protein_0672", "protein_0806", "protein_0807", "protein_0809", "protein_0810", "protein_0833", "protein_0932", "protein_0910")
dprot = load_observed_proteomics()

sum(filter(PR, id%in%rib_prot)$sim_mass)/sum(PR$sim_mass)
sum(filter(PR, id%in%rib_prot)$obs_mass)/sum(PR$obs_mass)

tail(X)

X = d_c$Protein/(d_c$RNA+d_c$Protein)
#plot(X, type="l")
d_c$RNA
