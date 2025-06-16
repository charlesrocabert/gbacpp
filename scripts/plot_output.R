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
  D              = d_b[i,-which(names(d_b)%in%c("condition", "h2o"))]
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
  X                = d_p[i, -which(names(d_p)%in%c("condition"))]
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
  reaction_ids   = names(d_p[,-which(names(d_p)%in%c("condition"))])
  excluded       = collect_excluded_proteins(reaction_ids, enz_comp)
  sim_proteomics = calculate_simulated_proteomics(model_path, model_name, d_p, i)
  Xsum           = sum(obs_proteomics$obs_mass)
  Xsim           = sum(filter(obs_proteomics, protein%in%rownames(sim_proteomics))$obs_mass)
  #print(paste("> Modeled proteome fraction", Xsim/Xsum))
  D              = data.frame(sim_proteomics$p_id, sim_proteomics$sim_mass, sim_proteomics$sim_fraction)
  names(D)       = c("id", "sim_mass", "sim_fraction")
  D              = filter(D, id%in%obs_proteomics$protein)
  D$obs_mass     = obs_proteomics[D$id, "obs_mass"]
  D              = D[!D$id%in%excluded,]
  D$obs_fraction = D$obs_mass/sum(D$obs_mass)
  D$sim_fraction = D$sim_mass/sum(D$sim_mass)#*0.36
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

### Build the growth law data ###
build_growth_law_data <- function( d_p )
{
  Phi = c()
  for(i in seq(1, dim(d_p)[1]))
  {
    Y    = d_p[i,-which(names(d_p)%in%c("condition"))]
    Ysum = rowSums(Y)
    Phi  = c(Phi, Y$Ribosome/Ysum)
  }
  D        = data.frame(d_p$condition, Phi)
  names(D) = c("condition", "Phi")
  return(D)
}

### Plot predicted protein fraction ###
plot_predicted_protein_fraction <- function( d_c, d_state )
{
  dl                  = d_c[,-which(names(d_c)%in%c("condition", "h2o"))]
  dl$sum              = rowSums(dl)
  dl$Protein_fraction = dl$Protein/dl$sum
  dl$condition        = d_c$condition
  dl$mu               = d_state$mu
  last_prot_fraction  = dl$Protein_fraction[dim(dl)[1]]
  p = ggplot(dl, aes(mu, Protein_fraction)) +
    geom_hline(yintercept=0.54727, color="pink") +
    geom_line() +
    scale_x_log10() +
    xlab("Glucose concentration (g/L)") +
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
    annotate("text", x=1e-5, y=5.5e-2, label=paste0("italic(r)^2", "==", round(r2,5)), hjust=0, parse=T) +
    xlab("Observed") +
    ylab("Simulated") +
    ggtitle("Metabolite mass fractions") +
    theme_classic()
  return(p)
}

### Plot predicted proteomics ###
plot_predicted_proteomics <- function( pr_data, R2, r2 )
{
  p = ggplot(pr_data, aes(obs_fraction, sim_fraction)) +
    geom_abline(slope=1, intercept=0, color="pink") +
    geom_line(aes(obs_fraction, obsp3), color="grey", lty=2) +
    geom_line(aes(obs_fraction, obsm3), color="grey", lty=2) +
    geom_line(aes(obs_fraction, obsp10), color="grey", lty=3) +
    geom_line(aes(obs_fraction, obsm10), color="grey", lty=3) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_x_log10() + scale_y_log10() +
    annotate("text", x=0.001, y=0.3, label=paste0("italic(R)^2", "==", round(R2,5)), hjust=0, parse=T) +
    annotate("text", x=0.001, y=0.1, label=paste0("italic(r)^2", "==", round(r2,5)), hjust=0, parse=T) +
    #geom_text_repel(aes(label=id), size = 3.5) +
    xlab("Observed") +
    ylab("Simulated") +
    ggtitle("Proteomics") +
    theme_classic()
  return(p)
}


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)

### Select the model ###
version     = "v1"
#output_path = "../output"
output_path = "/Users/charlesrocabert/Desktop/mmsyn_output"
model_path  = "../csv_models"
model_name  = paste0("mmsyn_fcr_", version)

### Load the experimental dataset ###
d_obs = read.table("/Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/data/wet_experiments/observed_mu.csv", sep=";", h=T)
Davg = tapply(d_obs$mu, d_obs$glc, mean)
Davg = data.frame(glc=as.numeric(names(Davg)), mu=Davg)

### Build the glucose vector ###
glc_vec = 10.0^seq(-5.9, 1, by=0.1)

### Load the model results ###
d_state     = read.table(paste0(output_path, "/", model_name, "_state_optimum.csv"), h=T, sep=";", check.names=F)
d_b         = read.table(paste0(output_path, "/", model_name,"_b_optimum.csv"), h=T, sep=";", check.names=F)
d_p         = read.table(paste0(output_path, "/", model_name,"_p_optimum.csv"), h=T, sep=";", check.names=F)
d_c         = read.table(paste0(output_path, "/", model_name,"_c_optimum.csv"), h=T, sep=";", check.names=F)
d_state$glc = glc_vec

### Build final datasets ###
MF_data  = build_mass_fractions_data(d_b, dim(d_b)[1])
PR_data  = build_proteomics_data(model_path, model_name, d_p, dim(d_p)[1])
MF_cor   = mass_fractions_correlation(d_b, dim(d_b)[1])
PR_cor   = proteomics_correlation(model_path, model_name, d_p, dim(d_p)[1])
MF_cor_D = data.frame()
PR_cor_D = data.frame()
for (i in seq(1, dim(d_b)[1]))
{
  res      = mass_fractions_correlation(d_b, i)
  MF_cor_D = rbind(MF_cor_D, res)
  res      = proteomics_correlation(model_path, model_name, d_p, i)
  PR_cor_D = rbind(PR_cor_D, res)
}
MF_cor_D$condition  = d_b$condition
MF_cor_D$glc        = glc_vec
names(MF_cor_D)     = c("r2", "pval", "R2", "condition", "glc")
PR_cor_D$condition  = d_p$condition
PR_cor_D$glc        = glc_vec
names(PR_cor_D)     = c("r2", "pval", "R2", "condition", "glc")
growth_law_data     = build_growth_law_data(d_p)
growth_law_data$glc = glc_vec
growth_law_data$mu  = d_state$mu

### Make plots ###
p1 = ggplot(d_state, aes(glc, mu)) +
  geom_line() +
  geom_point(data=Davg, aes(glc, mu, color="Observed")) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Glucose concentration (g/L)") +
  ylab("Growth rate") +
  ggtitle(paste0("Growth rate (max=",round(max(d_state$mu),3),")")) +
  theme_classic() +
  theme(legend.position="none")

p2 = plot_predicted_protein_fraction(d_c, d_state)

p3 = plot_predicted_mass_fractions(MF_data, MF_cor[3], MF_cor[1])

p4 = plot_predicted_proteomics(PR_data, PR_cor[3], PR_cor[1])

p5 = ggplot(MF_cor_D, aes(glc, R2)) +
  geom_line() +
  scale_x_log10() +
  xlab("Glucose concentration (g/L)") +
  ylab(TeX("$R^2$")) +
  ggtitle(paste0("Metabolite mass fraction correlation (max=",round(MF_cor[3],3),")")) +
  theme_classic()

rib_prot = c("protein_0025", "protein_0027", "protein_0082", "protein_0137", "protein_0148", "protein_0149", "protein_0198", "protein_0199", "protein_0238", "protein_0294", "protein_0362", "protein_0365", "protein_0422", "protein_0482", "protein_0499", "protein_0500", "protein_0501", "protein_0930", "protein_0526", "protein_0540", "protein_0637", "protein_0638", "protein_0644", "protein_0646", "protein_0647", "protein_0648", "protein_0653", "protein_0654", "protein_0655", "protein_0656", "protein_0657", "protein_0658", "protein_0659", "protein_0660", "protein_0661", "protein_0662", "protein_0663", "protein_0664", "protein_0665", "protein_0666", "protein_0667", "protein_0668", "protein_0669", "protein_0670", "protein_0671", "protein_0672", "protein_0806", "protein_0807", "protein_0809", "protein_0810", "protein_0833", "protein_0932", "protein_0910")
dprot    = load_observed_proteomics()
obs_phi  = sum(filter(dprot, protein%in%rib_prot)$obs_mass)/sum(dprot$obs_mass)
obs_phi
p6 = ggplot(growth_law_data, aes(mu, Phi)) + #*0.3627)) +
  geom_line() +
  geom_point(aes(x=0.34, y=obs_phi, color="Observed Phi")) +
  #geom_point(aes(x=0.4, y=obs_phi+0.16274/0.54727, color="Observed Phi + rRNA")) +
  xlab("Growth rate") +
  ylab(TeX("$\\Phi$ (weighted by 0.36)")) +
  ggtitle("Growth law") +
  #scale_y_log10() +
  theme_classic() +
  theme(legend.position="none")
p6

plot_grid(p1, p2, p3, p4, p5, p6, ncol=3, labels="AUTO")
# ggplot(PR_cor_D, aes(glc, R2)) +
#   scale_x_log10() +
#   geom_point()
