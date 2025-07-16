#################### Title:

# 3_merge_results.R
# Project: "SelectSim_Reviews"

#################### Author:

# Miljan Petrovic, postdoctoral researcher,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# March, 2025.

#################### Description:

# R script to merge results of different PM enrichment elastic net runs

#################### 
rm(list=ls())
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews/scripts/treatment_analysis")
library(ggplot2)
library(grid)
library(gridExtra)
slog <- function(y) sign(y)*log(1+abs(y))
params <- list()
params$cols = c("#66c2a5","#fc8d62","#8da0cb","#d9d57d","#4c6e1b","#7030a0","#946b2d","#0392cf","#f7cac9","#e78ac3","#a6d854","#961203","#AFE1AF","#DC123C","#00468b","#CF9FFF","#BFC0D7")
params$fsz = 16
save_ggplot_custom <- function(plot_name, plot_final, width=40, height=24) {
  ggsave(paste0("", plot_name, ".png"), plot=plot_final, device="png", width=width, height=height, units = "cm")
  return(NA)
}
save_ggplot_custom_pdf <- function(plot_name, plot_final, width=40, height=24) {
  ggsave(paste0("", plot_name, ".pdf"), plot=plot_final, device="pdf", width=width, height=height, units = "cm")
  return(NA)
}

# Load results
# load("../../results/treatment_analysis/EN_binomial_1.RData")
load("../../results/treatment_analysis/EN_binomial_co_1.RData")
# load("../../results/treatment_analysis/EN_binomial_2.RData")
load("../../results/treatment_analysis/EN_binomial_co_2.RData")
# load("../../results/treatment_analysis/EN_gaussian.RData")
load("../../results/treatment_analysis/EN_gaussian_co.RData")
# load("../../results/treatment_analysis/EN_binomial_1_drug.RData")
load("../../results/treatment_analysis/EN_binomial_co_1_drug.RData")
# load("../../results/treatment_analysis/EN_binomial_2_drug.RData")
load("../../results/treatment_analysis/EN_binomial_co_2_drug.RData")
# load("../../results/treatment_analysis/EN_gaussian_drug.RData")
load("../../results/treatment_analysis/EN_gaussian_co_drug.RData")

# Structure results (WITH DRUGS) in data tables
results <- EN_binomial_co_1_drug$coeff.matrix[ EN_binomial_co_1_drug$summary$`times selected`==10, ]
type_predictor <- rep("cancer_type", nrow(results))
type_predictor[ sapply(rownames(results), function(x) return( toupper(x)==x )) ] <- "gene"
type_predictor[ sapply(rownames(results), function(x) return( grepl(" - ",x))) ] <- "co-mutation"
type_predictor[118:nrow(results)] <- "treatment"
results_df_evidence <- data.frame( predictor = rep(rownames(results), 10), 
                          type_predictor = rep(type_predictor, 10),
                          coef = as.vector(as.matrix(results)),
                          repetition=seq(1,10)
)
results_df_evidence$analysis <- "binomial/met_status"

results <- EN_gaussian_co_drug$coeff.matrix[ EN_gaussian_co_drug$summary$`times selected`==10, ]
type_predictor <- rep("cancer_type", nrow(results))
type_predictor[ sapply(rownames(results), function(x) return( toupper(x)==x )) ] <- "gene"
type_predictor[ sapply(rownames(results), function(x) return( grepl(" - ",x))) ] <- "co-mutation"
type_predictor[171:nrow(results)] <- "treatment"
results_df_count <- data.frame( predictor = rep(rownames(results), 10), 
                                type_predictor = rep(type_predictor, 10),
                                coef = as.vector(as.matrix(results)),
                                repetition=seq(1,10)
)
results_df_count$analysis <- "gaussian/met_count"

results <- EN_binomial_co_2_drug$coeff.matrix[ EN_binomial_co_2_drug$summary$`times selected`==10, ]
type_predictor <- rep("cancer_type", nrow(results))
type_predictor[ sapply(rownames(results), function(x) return( toupper(x)==x )) ] <- "gene"
type_predictor[ sapply(rownames(results), function(x) return( grepl(" - ",x))) ] <- "co-mutation"
type_predictor[127:nrow(results)] <- "treatment"
results_df_nonzero_count <- data.frame( predictor = rep(rownames(results), 10), 
                                   type_predictor = rep(type_predictor, 10),
                                   coef = as.vector(as.matrix(results)),
                                   repetition=seq(1,10)
)
results_df_nonzero_count$analysis <- "binomial/met_count>0"

results_df_with_drugs <- rbind(results_df_evidence, results_df_count, results_df_nonzero_count)
results_df_with_drugs$treatment_included <- "yes"

# Structure results (WITHOUT DRUGS) in data tables
results <- EN_binomial_co_1$coeff.matrix[ EN_binomial_co_1$summary$`times selected`==10, ]
type_predictor <- rep("cancer_type", nrow(results))
type_predictor[ sapply(rownames(results), function(x) return( toupper(x)==x )) ] <- "gene"
type_predictor[ sapply(rownames(results), function(x) return( grepl(" - ",x))) ] <- "co-mutation"
results_df_evidence <- data.frame( predictor = rep(rownames(results), 10), 
                                   type_predictor = rep(type_predictor, 10),
                                   coef = as.vector(as.matrix(results)),
                                   repetition=seq(1,10)
)
results_df_evidence$analysis <- "binomial/met_status"

results <- EN_gaussian_co$coeff.matrix[ EN_gaussian_co$summary$`times selected`==10, ]
type_predictor <- rep("cancer_type", nrow(results))
type_predictor[ sapply(rownames(results), function(x) return( toupper(x)==x )) ] <- "gene"
type_predictor[ sapply(rownames(results), function(x) return( grepl(" - ",x))) ] <- "co-mutation"
results_df_count <- data.frame( predictor = rep(rownames(results), 10), 
                                type_predictor = rep(type_predictor, 10),
                                coef = as.vector(as.matrix(results)),
                                repetition=seq(1,10)
)
results_df_count$analysis <- "gaussian/met_count"

results <- EN_binomial_co_2$coeff.matrix[ EN_binomial_co_2$summary$`times selected`==10, ]
type_predictor <- rep("cancer_type", nrow(results))
type_predictor[ sapply(rownames(results), function(x) return( toupper(x)==x )) ] <- "gene"
type_predictor[ sapply(rownames(results), function(x) return( grepl(" - ",x))) ] <- "co-mutation"
results_df_nonzero_count <- data.frame( predictor = rep(rownames(results), 10), 
                                        type_predictor = rep(type_predictor, 10),
                                        coef = as.vector(as.matrix(results)),
                                        repetition=seq(1,10)
)
results_df_nonzero_count$analysis <- "binomial/met_count>0"

results_df_without_drugs <- rbind(results_df_evidence, results_df_count, results_df_nonzero_count)
results_df_without_drugs$treatment_included <- "no"

results_df <- rbind(results_df_with_drugs, results_df_without_drugs)

save(results_df, file="../../results/treatment_analysis/EN_results_summary.RData")


# ADD RESULTS INTO SUPPLEMENTARY TABLE 9
wb <- openxlsx::loadWorkbook("../../results/treatment_analysis/supplementary_table_9.xlsx")
openxlsx::addWorksheet(wb, "coefficients")
openxlsx::writeData( wb, 
                     sheet = "coefficients",
                     subset(results_df, analysis %in% c("binomial/met_count>0","gaussian/met_count") & treatment_included=="yes" ) )
openxlsx::saveWorkbook(wb, "../../results/treatment_analysis/supplementary_table_9.xlsx", overwrite = TRUE)


