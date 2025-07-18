#################### Title:

# normal_mutations/1_preprocess_Abby2023.R
# Project: "SelectSim_Reviews"

#################### Author:

# Miljan Petrovic, postdoctoral researcher,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# April, 2025.

#################### Description:

# R script to preprocess skin mutation data from Fowler2021 
rm(list=ls())
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews")
library(Matrix)
source("scripts/rerun_selectSim_new/routines_coincid_test.R")

# mutset <- "allmut"
mutset <- "cancermut"
#################### 

##### ABBY ESOPHAGUS TISSUE #####

# Load mutations from normal tissue
# normal_samples <- c("PD37449","PD34201","PD34200","PD37590","PD37266","PD34199","PD28690") 
SuppT <- readxl::read_xlsx("data/raw/normal_mutations/Abby2023/41588_2022_1280_MOESM4_ESM.xlsx", sheet=3, skip=27, col_names = T)
SuppT <- SuppT[1:630,]

SuppT$subject <- sapply(SuppT$Sample, function(x) stringr::str_extract(x, "^[0-9A-Z]*") )
SuppT$sample <- SuppT$Sample
SuppT$gene <- SuppT$Gene
input_maf <- SuppT

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
  'ignore' = c("silent"),
  'truncating'= c('splice_region','nonsense','inframe','frameshift','ess_splice'),
  'missense' = c('missense')
)

input_maf$summary <- input_maf$Effect
saveRDS(input_maf, file="data/preprocessed/normal_mutations/Abby2023_maf.rds")

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
maf_m$mut_filter <- paste0( maf_m$gene, ":", sapply(maf_m$Protein, function(x) substr(x,3,nchar(x)-1) ) )
if (mutset == "cancermut") variant_catalogue$mut_filter <- paste0( variant_catalogue$gene, ":", variant_catalogue$mut )
if (mutset == "cancermut") maf_m <- subset( maf_m, mut_filter %in% subset(variant_catalogue, oncogenic %in% c("Oncogenic","Likely Oncogenic"))$mut_filter )
GAM_m <- Matrix::sparseMatrix( i = match(maf_m$gene, genes),
                               j = match(maf_m$sample, samples_all),
                               x = rep(1, nrow(maf_m)), dims = c(length(genes), length(samples_all)),
                               dimnames = list(genes, samples_all) )
GAM_m_bin <- GAM_m; GAM_m_bin@x[GAM_m_bin@x<1] <- 0; GAM_m_bin@x[GAM_m_bin@x>=1] <- 1

if (mutset == "allmut") {
  genes_ind <- which(Matrix::rowSums(GAM_t_bin) + Matrix::rowSums(GAM_m_bin) > 5)
} else {
  genes_ind <- which(Matrix::rowSums(GAM_t_bin) + Matrix::rowSums(GAM_m_bin) > 0)
}
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

# Calculate template matrix
# temp_mat_t <- matrix(0, nrow(GAM_t_bin), 0)
# temp_mat_m <- matrix(0, nrow(GAM_m_bin), 0)
# for (iter_block in 1:length(blocks)) {
#   inds <- blocks[[iter_block]]
#   
#   gfreq_t <- as.vector(Matrix::rowSums(GAM_t_bin[,inds]) / length(inds))
#   sfreq_t <- as.vector(length(inds) * as.numeric(tmb_t[inds]) / sum(as.numeric(tmb_t[inds])))
#   block_temp <- outer( gfreq_t, sfreq_t ); block_temp[block_temp>1] <- 1; dimnames(block_temp) <- list(rownames(GAM_m_bin),inds)
#   temp_mat_t <- cbind(temp_mat_t, block_temp)
#   
#   gfreq_m <- as.vector(Matrix::rowSums(GAM_m_bin[,inds]) / length(inds))
#   sfreq_m <- as.vector(length(inds) * as.numeric(tmb_m[inds]) / sum(as.numeric(tmb_m[inds])))
#   block_temp <- outer( gfreq_m, sfreq_m ); block_temp[block_temp>1] <- 1; dimnames(block_temp) <- list(rownames(GAM_m_bin),inds)
#   temp_mat_m <- cbind(temp_mat_m, block_temp)
# }
# temp_mat_t <- temp_mat_t[,colnames(GAM_t_bin)]
# temp_mat_m <- temp_mat_m[,colnames(GAM_m_bin)]
# temp_mat_total <- pmax(temp_mat_t, temp_mat_m)
# 
# # Get sample weights
# block_means <- lapply(blocks, function(x) return( median(tmb_t[x]+tmb_m[x])) )
# FCs <- (tmb_t+tmb_m) / mean(unlist(block_means)); FCs[FCs<=1] <- 1
# sample_weights <- weighting(FCs, lam = 0.3, tau = 1)
# 
# # Build nulls
# set.seed(42)
# nulls <- list()
# nulls <- lapply(seq(1,1000), function(x) {
#   res <- gen_null_init( temp_mat=temp_mat_total, 
#                         observed_mat=GAM_mut_bin, 
#                         type="rejection_row_correction", 
#                         blocks=lapply(blocks, function(x) return(match(x, colnames(GAM_mut_bin))) )
#   )
#   return( res )
# })
# outliers <- retrieveOutliers(nulls, Matrix::rowSums(GAM_mut_bin), Matrix::colSums(GAM_mut_bin)); nulls <- nulls[!outliers]

# Save preprocessed
# run_data <- list()
# run_data$GAM_t_bin <- GAM_t_bin
# run_data$GAM_m_bin <- GAM_m_bin
# run_data$GAM_mut_bin <- GAM_mut_bin
# run_data$blocks <- blocks
# run_data$temp_mat_total <- temp_mat_total
# run_data$sample_weights <- sample_weights
# run_data$tmb_t <- tmb_t
# run_data$tmb_m <- tmb_m
# run_data$nulls <- nulls
# run_data$blocks <- blocks
# run_data$sample2patient <- sample2patient
# run_data$samples_all <- samples_all
# save(run_data, file="results/normal_mutations/run_data_Fowler2021.RData")

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
save(run_data_selectSim, file=paste0("results/normal_mutations/run_data_Abby2023_forSelectSim_",mutset,".RData"))






