# Author: Miljan Petrovic
rm(list = ls())
#setwd("/Users/mpetrov2/Documents/POSTDOC/UNIL/ArvindProject/BladderAnalysis_reprod")
library(SelectX)
library(dplyr)
library(tictoc)
library(ggplot2)
library(gridExtra)

##################################################
# LOAD DATA
##################################################

tcga_rundata_full <- readRDS(file="./data/tcga_BLCA.rds")
seeds <- c(95, 90, 37, 13, 40, 56, 32, 45, 69, 6)
li_cancer_rundata <- readRDS(file="./data/li_et_al_cancer_run_data_v3.rds")
li_normal_rundata <- readRDS(file="./data/li_et_al_normal_run_data_v3.rds")
li_cancer_rundata$M$tmb$missense$mutation <- as.numeric( li_cancer_rundata$M$tmb$missense$mutation )
li_cancer_rundata$M$tmb$truncating$mutation <- as.numeric( li_cancer_rundata$M$tmb$truncating$mutation )
li_normal_rundata$M$tmb$missense$mutation <- as.numeric( li_normal_rundata$M$tmb$missense$mutation )
li_normal_rundata$M$tmb$truncating$mutation <- as.numeric( li_normal_rundata$M$tmb$truncating$mutation )

##################################################
# RUN SELECTSIM on TCGA BLADDER (MULTIPLE TIMES)
##################################################
result_obj_tcga <- as.list(rep(NA,length(seeds)))
tic()
for (iter in 1:length(seeds)) {
  
  # subsample tcga run data
  set.seed(seeds[iter]); inds_tcga <- sample(398,113)
  tcga_rundata <- tcga_rundata_full
  tcga_rundata$M$M$missense <- tcga_rundata_full$M$M$missense[,inds_tcga]
  tcga_rundata$M$M$truncating <- tcga_rundata_full$M$M$truncating[,inds_tcga]
  tcga_rundata$M$tmb$missense <- tcga_rundata_full$M$tmb$missense[inds_tcga,]
  tcga_rundata$M$tmb$truncating <- tcga_rundata_full$M$tmb$truncating[inds_tcga,]
  tcga_rundata$sample.class <- tcga_rundata_full$sample.class[inds_tcga]
  
  # run selectSim
  result_obj_tcga[[iter]] <- selectX(  M = tcga_rundata$M,
                                       sample.class = tcga_rundata$sample.class,
                                       alteration.class = tcga_rundata$alteration.class,
                                       n.cores = 1,
                                       min.freq = ceiling(0.001*length(tcga_rundata$sample.class)),
                                       n.permut = 1000,
                                       lambda = 0.3,
                                       tao = 1,
                                       save.object = FALSE,
                                       verbose = FALSE,
                                       estimate_pairwise = FALSE,
                                       maxFDR = 0.25,
                                       seed = seeds[iter])
  result_obj_tcga[[iter]]$result$sample_inds <- rep(paste0(inds_tcga,collapse = ","), times=nrow(result_obj_tcga[[iter]]$result))
}
saveRDS(result_obj_tcga, file = "./results/tcga_results_obj_with_nulls_bladder_merged.rds")
toc()
tcga_results_bladder <- lapply(result_obj_tcga, function(x) x$result)
nrows_tcga <- unlist(lapply(tcga_results_bladder, function(x) nrow(x)))
tcga_results_bladder_merged <- do.call("rbind", tcga_results_bladder)
tcga_results_bladder_merged$seed <- rep(seeds, times=nrows_tcga)
saveRDS(tcga_results_bladder_merged, file = "./results/tcga_results_bladder_merged.rds")


##################################################
# RUN SELECTSIM on Li BLADDER
##################################################

result_obj_li_cancer <- selectX(  M = li_cancer_rundata$M,
                                     sample.class = li_cancer_rundata$sample.class,
                                     alteration.class = li_cancer_rundata$alteration.class,
                                     n.cores = 1,
                                     min.freq = ceiling(0.001*length(li_cancer_rundata$sample.class)),
                                     n.permut = 1000,
                                     lambda = 0.3,
                                     tao = 1,
                                     save.object = FALSE,
                                     verbose = FALSE,
                                     estimate_pairwise = FALSE,
                                     maxFDR = 0.25,
                                     seed = 42)
result_obj_li_normal <- selectX(  M = li_normal_rundata$M,
                                  sample.class = li_normal_rundata$sample.class,
                                  alteration.class = li_normal_rundata$alteration.class,
                                  n.cores = 1,
                                  min.freq = ceiling(0.001*length(li_normal_rundata$sample.class)),
                                  n.permut = 1000,
                                  lambda = 0.3,
                                  tao = 1,
                                  save.object = FALSE,
                                  verbose = FALSE,
                                  estimate_pairwise = FALSE,
                                  maxFDR = 0.25,
                                  seed = 42)
result_obj_li_normal_noTMB <- selectX(  M = li_normal_rundata$M,
                                        sample.class = li_normal_rundata$sample.class,
                                        alteration.class = li_normal_rundata$alteration.class,
                                        n.cores = 1,
                                        min.freq = ceiling(0.001*length(li_normal_rundata$sample.class)),
                                        n.permut = 1000,
                                        lambda = 0, #weighting function of TMB fixed to a constant 1
                                        tao = 1,
                                        save.object = FALSE,
                                        verbose = FALSE,
                                        estimate_pairwise = FALSE,
                                        maxFDR = 0.25,
                                        seed = 42)
#sanity check on using TMB for normal samples
ff = subset(result_obj_li_normal$result, FDR==TRUE)
ff_noTMB = subset(result_obj_li_normal_noTMB$result, FDR==TRUE)
table( ff$overlap )
table( ff_noTMB$overlap )

li_results_bladder_merged <- rbind(result_obj_li_cancer$result, result_obj_li_normal$result, result_obj_li_normal_noTMB$result)
li_results_bladder_merged$cohort <- rep(c("Cancer","Normal","Normal"), times=c(nrow(result_obj_li_cancer$result),nrow(result_obj_li_normal$result),nrow(result_obj_li_normal_noTMB$result)))
li_results_bladder_merged$TMB_accounted <- rep(c("Yes","Yes","No"), times=c(nrow(result_obj_li_cancer$result),nrow(result_obj_li_normal$result),nrow(result_obj_li_normal_noTMB$result)))
save(result_obj_li_cancer, result_obj_li_normal, result_obj_li_normal_noTMB, file = "./results/li_obj_bladder_merged.RData")
saveRDS(li_results_bladder_merged, file = "./results/li_results_bladder_merged.rds")