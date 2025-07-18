#################### Title:

# normal_mutations/1_preprocess_Fowler.R
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

#################### 

##### LAWSON_2020 UROTHELIUM #####

# Load mutations from normal tissue
input_maf <- readxl::read_xlsx("data/raw/normal_mutations/Lawson2020/aba8347_tables3.xlsx", sheet=1, skip=1, col_names = T) #targeted sequencing
input_maf <- subset(input_maf, histological_feature %in% c("Urothelium","Tumour","CIS") )

# Build vcf file to call variants online
# vcf <- cbind( input_maf[,c("chr","pos")], data.frame("ID"=rep(".",nrow(input_maf))), input_maf[,c("ref","mut")], data.frame("QUAL"=rep(30,nrow(input_maf)), "FILTER"=rep("PASS",nrow(input_maf))) )
# colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER")
# AF <- paste0("AF=", input_maf$vaf)
# NS <- paste0("NS=", input_maf$mut_depth)
# DP <- paste0("DP=", input_maf$total_depth)
# vcf$INFO <- paste0(AF,";",NS,";",DP)
# vcf$SAMPLE <- input_maf$microbiopsy_id
# vcf <- vcf[order(vcf$POS),]
# write.table(vcf, file="data/preprocessed/normal_mutations/lawson.vcf", sep="\t", quote=F, row.names=F, col.names=T) # for: https://run.opencravat.org/ on hg19 + ClinVar

# Load vcf file transfer annotations of variant consequences
post_vcf <- read.delim("data/preprocessed/normal_mutations/lawson_vep_output.txt",skip=0)
post_vcf <- subset(post_vcf, Gene !="-")
mutation_type = list(
  'truncating'= c('stop_gained','inframe_deletion','inframe_deletion,NMD_transcript_variant','stop_gained,splice_region_variant','stop_gained,NMD_transcript_variant','stop_gained,inframe_deletion','frameshift_variant','start_lost','frameshift_variant,splice_region_variant','frameshift_variant,NMD_transcript_variant','stop_gained,inframe_deletion'),
  'missense' = c('missense_variant','missense_variant,NMD_transcript_variant','missense_variant,splice_region_variant','missense_variant,splice_region_variant,NMD_transcript_variant')
)
post_vcf <- subset(post_vcf, Consequence %in% unlist(mutation_type))
# location2consequence <- split(post_vcf$Consequence, post_vcf$Location); location2consequence <- lapply(location2consequence, unique)
# ii <- which(unlist(lapply(location2consequence,length))==2)
# xx <- location2consequence[ii]
#location2consequence[ii] <- "stop_gained"
# post_vcf <- subset(post_vcf, !(Location %in% names(xx) & Consequence!="stop_gained"))
# location2consequence <- split(post_vcf$Consequence, post_vcf$Location); location2consequence <- lapply(location2consequence, unique)

# Assign consequence of each variant (mapping variants by location)
input_maf$Location <- paste0(input_maf$chr,":",input_maf$pos,"-", input_maf$pos)
maf <- post_vcf
maf <- maf[!duplicated(maf),]
maf$microbiopsy_id <- NA
maf$donor <- NA
maf$histology <- NA
maf$vaf <- NA
format_string <- function(input) {
  matches <- regmatches(input, regexec("^([A-Za-z0-9]+):(\\d+)-\\d+$", input))[[1]]
  if (length(matches) == 0) return(NA)
  paste0(matches[2], ":", matches[3], "-", matches[3])
}
maf$Location_dummy <- sapply(maf$Location, format_string)
all_locs <- unique(maf$Location_dummy)
for (iter_loc in 1:length(all_locs)) {
  loc <- all_locs[iter_loc]
  tmp <- subset(input_maf, Location == loc)
  inds <- which(maf$Location_dummy == loc)
  maf$samples[inds] <- paste0( tmp$microbiopsy_id, collapse=";" )
  maf$vaf[inds] <- paste0( tmp$vaf, collapse=";" )
  maf$histology[inds] <- paste0( tmp$histological_feature, collapse=";" )
}
maf <- maf[,c("Location","Location_dummy","SYMBOL","Consequence","Amino_acids","Protein_position","samples","vaf","histology")]
colnames(maf) <- c("Location","Location_dummy","gene","Consequence","Amino_acids","Protein_position","samples","vaf","histology")
maf <- maf[!duplicated(maf),]
maf <- maf[!duplicated(maf[,c("Location_dummy","gene","Consequence","samples")]),]
maf$mut_filter <- paste0( maf$gene, ":", sapply(maf$Amino_acids, function(x) substr(x,1,1)), maf$Protein_position  )

# Unravel maf which now has one row per mutation into one row per pair of mutation and sample
maf_final <- data.frame()
for (iter in 1:nrow(maf)) {
  samples <- strsplit(maf$samples[iter], ";")[[1]]
  histologies <- strsplit(maf$histology[iter], ";")[[1]]
  vafs <- strsplit(maf$vaf[iter], ";")[[1]]
  tmp_df <- maf[rep(iter, length(samples)),]
  tmp_df$samples <- samples
  tmp_df$histology <- histologies
  tmp_df$vaf <- vafs
  maf_final <- rbind(maf_final, tmp_df)
}
maf <- maf_final
colnames(maf) <- c("Location","Location_dummy","gene","Consequence","Amino_acids","Protein_position","sample","vaf","histology","mut_filter")
maf$subject <- sapply(maf$sample, function(x) paste0( strsplit(x, "_")[[1]][1:2], collapse="_") )
maf$subject_type <- "non-cancer"
maf$subject_type[maf$subject %in% c("C01_49M","C04_72M","C05_75M")] <- "cancer-non-muscle-invasive"
maf$subject_type[maf$subject %in% c("C02_61F","C03_67M")] <- "cancer-muscle-invasive"
maf$AA_change <- paste0(maf$gene, ":", sapply(maf$Amino_acids, function(x) substr(x,1,1)), maf$Protein_position, sapply(maf$Amino_acids, function(x) substr(x,3,3)) )
saveRDS(maf, file="data/preprocessed/normal_mutations/lawson2020_maf.rds")


# BUILD SELECT_SIM RUN DATA
# mutset <- "allmut"
mutset <- "cancermut"
maf_full <- readRDS("data/preprocessed/normal_mutations/lawson2020_maf.rds")
mutation_type = list(
  'truncating'= c('stop_gained','inframe_deletion','inframe_deletion,NMD_transcript_variant','stop_gained,splice_region_variant','stop_gained,NMD_transcript_variant','stop_gained,inframe_deletion','frameshift_variant','start_lost','frameshift_variant,splice_region_variant','frameshift_variant,NMD_transcript_variant','stop_gained,inframe_deletion'),
  'missense' = c('missense_variant','missense_variant,NMD_transcript_variant','missense_variant,splice_region_variant','missense_variant,splice_region_variant,NMD_transcript_variant')
)
runs <- c("normal_samples_all_donors", "normal_samples_normal_donors", "normal_samples_cancer_donors","cancer_samples_cancer_donors")


run_data <- list()
for (iter_run in 1:length(runs)) {
  run <- runs[iter_run]
  if (run == "normal_samples_all_donors") {
    maf <- subset(maf_full, histology =="Urothelium")
  } else if (run == "normal_samples_normal_donors") {
    maf <- subset(maf_full, histology =="Urothelium" & subject_type == "non-cancer")
  } else if (run == "normal_samples_cancer_donors") {
    maf <- subset(maf_full, histology =="Urothelium" & subject_type != "non-cancer")
  } else if (run == "cancer_samples_cancer_donors") {
    maf <- subset(maf_full, histology !="Urothelium" & subject_type != "non-cancer")
  }
  
  patient2sample <- split(maf$sample, maf$subject); patient2sample <- lapply(patient2sample, unique)
  sample2patient <- split(maf$subject, maf$sample); sample2patient <- lapply(sample2patient, unique); sample2patient <- unlist(sample2patient[sort(names(sample2patient))])
  samples_all <- names(sample2patient)
  if (mutset == "cancermut") {data(oncokb_genes, package = "SelectSim"); oncokb_genes =  sort(oncokb_genes)}
  if (mutset == "cancermut") {data(variant_catalogue, package = "SelectSim"); variant_catalogue = variant_catalogue}
  if (mutset == "cancermut") {genes <- sort(intersect(unique(maf$gene), oncokb_genes))}
  if (mutset == "allmut") genes <- sort(unique(maf$gene))
  blocks <- split(names(sample2patient), sample2patient)
  
  # build truncating
  maf_t <- subset( maf, Consequence %in% mutation_type$truncating )
  maf_t <- subset( maf_t, gene %in% genes ) 
  GAM_t <- Matrix::sparseMatrix( i = match(maf_t$gene, genes),
                                 j = match(maf_t$sample, samples_all),
                                 x = rep(1, nrow(maf_t)), dims = c(length(genes), length(samples_all)),
                                 dimnames = list(genes, samples_all) )
  GAM_t_bin <- GAM_t; GAM_t_bin@x[GAM_t_bin@x<1] <- 0; GAM_t_bin@x[GAM_t_bin@x>=1] <- 1
  
  # build missense
  maf_m <- subset( maf, Consequence %in% mutation_type$missense )
  maf_m <- subset( maf_m, gene %in% genes )
  if (mutset == "cancermut") variant_catalogue$mut_filter <- paste0( variant_catalogue$gene, ":", variant_catalogue$mut )
  if (mutset == "cancermut") maf_m <- subset( maf_m, mut_filter %in% subset(variant_catalogue, oncogenic %in% c("Oncogenic","Likely Oncogenic"))$mut_filter)
  GAM_m <- Matrix::sparseMatrix( i = match(maf_m$gene, genes),
                                 j = match(maf_m$sample, samples_all),
                                 x = rep(1, nrow(maf_m)), dims = c(length(genes), length(samples_all)),
                                 dimnames = list(genes, samples_all) )
  GAM_m_bin <- GAM_m; GAM_m_bin@x[GAM_m_bin@x<1] <- 0; GAM_m_bin@x[GAM_m_bin@x>=1] <- 1
  
  GAM_full_Lawson <- pmax(GAM_m_bin,GAM_t_bin)
  save(GAM_full_Lawson, file=paste0("results/normal_mutations/GAM_Lawson2020_",run,".RData"))
  
  if (mutset == "cancermut") genes_ind <- which(Matrix::rowSums(GAM_t_bin) + Matrix::rowSums(GAM_m_bin) >= 2)
  if (mutset == "allmut") genes_ind <- which(Matrix::rowSums(GAM_t_bin) + Matrix::rowSums(GAM_m_bin) >= 5)
  GAM_t_bin <- GAM_t_bin[genes_ind,]
  GAM_m_bin <- GAM_m_bin[genes_ind,]
  genes <- genes[genes_ind]
  GAM_mut_bin <- GAM_t_bin+GAM_m_bin; GAM_mut_bin@x[GAM_mut_bin@x<1] <- 0; GAM_mut_bin@x[GAM_mut_bin@x>=1] <- 1
  
  tt <- subset( maf[!duplicated(maf),], Consequence %in% mutation_type$truncating )
  tmb_t <- rep(0, length(samples_all)); names(tmb_t) <- samples_all
  table_tmp <- table(tt$sample)
  tmb_t[names(table_tmp)] <- as.numeric(table_tmp); if (any(is.na(tmb_t))) tmb_t[is.na(tmb_t)] <- 0
  mm <- subset( maf[!duplicated(maf),], Consequence %in% mutation_type$missense )
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
  # save(run_data_selectSim, file=paste0("results/normal_mutations/run_data_Lawson2020_",run,".RData"))
  
    
}

























