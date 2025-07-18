#################### Title:

# normal_mutations/1_preprocess_Martincorena2018.R
# Project: "SelectSim_Reviews"

#################### Author:

# Miljan Petrovic, postdoctoral researcher,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# April, 2025.

#################### Description:

# R script to preprocess skin mutation data from Martincorena2018.R
# https://www.science.org/doi/10.1126/science.aau3879 

rm(list=ls())
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews")
library(Matrix)
# source("scripts/rerun_selectSim_new/routines_coincid_test.R")

# mutset <- "allmut"
mutset <- "cancermut"
#################### 

##### MARTINCORENA_2018 SKIN TISSUE #####

# Load mutations from normal tissue
input_maf <- readxl::read_xlsx("data/raw/normal_mutations/Martincorena2018/aau3879_tables2.xlsx", sheet=1, skip=16, col_names = T)
colnames(input_maf)[1] <- "sample"
input_maf$subject <- sapply(input_maf$sample, function(x) return(substr(x,1,7)) ) 
colnames(input_maf)[9] <- "summary"

# Number of samples/patients total
patient2sample <- split(input_maf$sample, input_maf$subject); patient2sample <- lapply(patient2sample, unique)
sample2patient <- split(input_maf$subject, input_maf$sample); sample2patient <- lapply(sample2patient, unique); sample2patient <- unlist(sample2patient[sort(names(sample2patient))])
samples_all <- names(sample2patient)
if (mutset == "cancermut") {data(oncokb_genes, package = "SelectSim"); oncokb_genes =  sort(oncokb_genes)}
if (mutset == "cancermut") {data(variant_catalogue, package = "SelectSim"); variant_catalogue = variant_catalogue}
if (mutset == "cancermut") genes <- sort(intersect(unique(input_maf$gene), oncokb_genes))
if (mutset == "allmut") genes <- sort(unique(input_maf$gene))
blocks <- split(names(sample2patient), sample2patient)

mutation_type = list(
  'ignore' = c("Synonymous","no-SNV"),
  'truncating'= c('Inframe Deletion','Stop_loss','Nonsense','Essential_Splice'),
  'missense' = c('Missense')
)

saveRDS(input_maf, file="data/preprocessed/normal_mutations/Martincorena2018_maf.rds")

# build truncating
maf_t <- subset( input_maf, summary %in% mutation_type$truncating )
maf_t <- subset( maf_t, gene %in% genes ) # "CMTR2"  "SETBP1" are not in "oncokb_truncating_genes" but have nonzeros in truncating gam of the package
GAM_t <- Matrix::sparseMatrix( i = match(maf_t$gene, genes),
                               j = match(maf_t$sample, samples_all),
                               x = rep(1, nrow(maf_t)), dims = c(length(genes), length(samples_all)),
                               dimnames = list(genes, samples_all) )
GAM_t_bin <- GAM_t; GAM_t_bin@x[GAM_t_bin@x<1] <- 0; GAM_t_bin@x[GAM_t_bin@x>=1] <- 1

# build missense
maf_m <- subset( input_maf, summary %in% mutation_type$missense )
maf_m <- subset( maf_m, gene %in% genes )
maf_m$mut_filter <- paste0( maf_m$gene, ":", sapply(maf_m$aachange, function(x) substr(x,1,nchar(x)-1) ) )
if (mutset == "cancermut") variant_catalogue$mut_filter <- paste0( variant_catalogue$gene, ":", variant_catalogue$mut )
if (mutset == "cancermut") maf_m <- subset( maf_m, mut_filter %in% subset(variant_catalogue, oncogenic %in% c("Oncogenic","Likely Oncogenic"))$mut_filter )
GAM_m <- Matrix::sparseMatrix( i = match(maf_m$gene, genes),
                               j = match(maf_m$sample, samples_all),
                               x = rep(1, nrow(maf_m)), dims = c(length(genes), length(samples_all)),
                               dimnames = list(genes, samples_all) )
GAM_m_bin <- GAM_m; GAM_m_bin@x[GAM_m_bin@x<1] <- 0; GAM_m_bin@x[GAM_m_bin@x>=1] <- 1

GAM_full_Martincorena18 <- pmax(GAM_m_bin,GAM_t_bin)
save(GAM_full_Martincorena18, file=paste0("results/normal_mutations/GAM_Martincorena2018_",mutset,".RData"))

genes_ind <- which(Matrix::rowSums(GAM_t_bin) + Matrix::rowSums(GAM_m_bin) >= 2)
GAM_t_bin <- GAM_t_bin[genes_ind,]
GAM_m_bin <- GAM_m_bin[genes_ind,]
genes <- genes[genes_ind]
GAM_mut_bin <- GAM_t_bin+GAM_m_bin; GAM_mut_bin@x[GAM_mut_bin@x<1] <- 0; GAM_mut_bin@x[GAM_mut_bin@x>=1] <- 1

# Calculate TMB values (original)
tt <- subset( input_maf[!duplicated(input_maf),], summary %in% mutation_type$truncating )
tmb_t <- rep(0, length(samples_all)); names(tmb_t) <- samples_all
table_tmp <- table(tt$sample)
tmb_t[names(table_tmp)] <- as.numeric(table_tmp); if (any(is.na(tmb_t))) tmb_t[is.na(tmb_t)] <- 0
mm <- subset( input_maf[!duplicated(input_maf),], summary %in% mutation_type$missense )
tmb_m <- rep(0, length(samples_all)); names(tmb_m) <- samples_all
table_tmp <- table(mm$sample)
tmb_m[names(table_tmp)] <- as.numeric(table_tmp); if (any(is.na(tmb_m))) tmb_m[is.na(tmb_m)] <- 0

# Restructure run data for selectSim package
run_data_selectSim <- list()
run_data_selectSim$M <- list()
run_data_selectSim$M$M <- list()
run_data_selectSim$M$M$missense <- as.matrix(GAM_m_bin)
run_data_selectSim$M$M$truncating <- as.matrix(GAM_t_bin)
run_data_selectSim$M$tmb <- list()
run_data_selectSim$M$tmb$missense <- data.frame( sample = samples_all, mutation = tmb_m, row.names = samples_all )
run_data_selectSim$M$tmb$truncating <- data.frame( sample = samples_all, mutation = tmb_t, row.names = samples_all )
run_data_selectSim$sample.class <- sample2patient
run_data_selectSim$alteration.class <- setNames(rep("MUT", length(genes)), genes)
save(run_data_selectSim, file=paste0("results/normal_mutations/run_data_Martincorena2018_forSelectSim_",mutset,".RData"))






