#################### Title:

# 1_preprocessTreatment.R
# Project: "SelectSim_Reviews"

#################### Author:

# Miljan Petrovic, postdoctoral researcher,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# March, 2025.

#################### Description:

# R script to preprocess drug treatment data of MSK 

#################### 
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews/scripts/treatment_analysis")

# RAW TREATMENT DATA DOWNLOADED FROM MSK-CHORD: https://www.cbioportal.org/study/summary?id=msk_chord_2024

# Load raw data and merge
drug_files <- list.files("../../data/raw/MSK_Treatment", full.names = TRUE)
for (iter in 1:length(drug_files)) {
  tmp <- strsplit(drug_files[iter], "/")[[1]]
  drug_name <- substr( tmp[length(tmp)], 1, nchar(tmp[length(tmp)])-4 )
  if (iter == 1) {
    drug_data <- read.delim(drug_files[iter], header=T)
    drug_data$drug <- drug_name
  } else {
    tmp_data <- read.delim(drug_files[iter], header=T)
    tmp_data$drug <- drug_name
    drug_data <- rbind(drug_data, tmp_data)
  }
  print(paste0(iter, "/", length(drug_files), " drugs loaded."))
}

# Save drug data
save(drug_data, file="../../data/raw/selectsim_analysis/analysis_data/primary_met_elastic_net/drug_data.RData")


