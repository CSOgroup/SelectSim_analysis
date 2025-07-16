#################### Title:

# 2_run_EN.R
# Project: "SelectSim_Reviews"

#################### Author:

# Debora Sesia, Miljan Petrovic,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# March, 2025.

#################### Description:

# R script to rerun PM enrichment elastic net together with drug treatment data of MSK 

#################### 
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews/scripts/treatment_analysis")
library(Matrix)
full2sparse <- function(mat) {
  inds <- which(mat!=0, arr.ind=T); inds <- cbind(inds, mat[inds])
  smat <- Matrix::sparseMatrix(i=inds[,1], 
                               j=inds[,2],
                               x=inds[,3],
                               dims = dim(mat), dimnames = list(rownames(mat),colnames(mat)),
                               symmetric = F
  )
  
  return(smat)
}
set.seed(42)

# Load data and code
# load("../../data/raw/selectsim_analysis/analysis_data/primary_met_elastic_net/PM_comutations(OLD).RData")
PM_comutations <- readRDS("../../data/raw/selectsim_analysis/analysis_data/primary_met_elastic_net/msk_p_m_fold_change_data_frame.rds")
PM_comutations <- subset(PM_comutations, category %in% c("M_specific","P_specific","PM_enriched","PM_switch"))
clinical_info=read.delim("../../data/raw/selectsim_analysis/analysis_data/primary_met_elastic_net/msk_full_clinical_table.txt")
load("../../data/raw/selectsim_analysis/analysis_data/primary_met_elastic_net/full_gam.RData")
load("../../data/raw/selectsim_analysis/analysis_data/primary_met_elastic_net/drug_data.RData")
source("EN_functions.R")

#-------------------------------------------------------------------------------

GAM=t(full_gam)
rm(full_gam)
dim(GAM)
rownames(GAM)=gsub("GENIE-MSK-","",rownames(GAM))

sample_to_keep=intersect(rownames(GAM),clinical_info$sample_id)
GAM=GAM[which(rownames(GAM)%in%sample_to_keep),]
clinical_info=clinical_info[which(clinical_info$sample_id%in%sample_to_keep),]

GAM=GAM[order(rownames(GAM)),]
clinical_info=clinical_info[order(clinical_info$sample_id),]

# add co-alteration variables (addition of interaction terms)
GAM_co=matrix(NA,ncol = nrow(PM_comutations),nrow = nrow(GAM))
colnames(GAM_co)=PM_comutations$id
rownames(GAM_co)=rownames(GAM)

temp=data.frame(apply(GAM_co, 2, function(x) sum(x)))

for(i in 1:nrow(PM_comutations)){
  co=PM_comutations$id[i]
  gene1=unlist(strsplit(co," - "))[1]
  gene2=unlist(strsplit(co," - "))[2]
  co_variable=GAM[,gene1]*GAM[,gene2]
  #length(which(co_variable==1))
  GAM_co[,i]=co_variable
} 

# add to GAM the tumor type 
GAM_tt=model.matrix( ~ 0 + cancer_type, clinical_info)
colnames(GAM_tt)=gsub("cancer_type","",colnames(GAM_tt))
rownames(GAM_tt)=clinical_info$sample_id

GAM_tt=GAM_tt[order(rownames(GAM_tt)),]


X=cbind(GAM_tt,GAM)
X_co=cbind(X,GAM_co)
X <- full2sparse(X)
X_co <- full2sparse(X_co)

# RUN PREDICTION WITHOUT DRUGS

#-------------------------------------------------------------------------------
# 1) binomial regression (logistic, y binary) with variable "is_evidence_of_mets"
#-------------------------------------------------------------------------------
y1=clinical_info$is_evidence_of_mets

EN_binomial_1=ENanalysis(y1,X,n.rep=10,type_regression="binomial")
save(EN_binomial_1,file="../../results/treatment_analysis/EN_binomial_1.RData")
#-------------------------------------------------------------------------------
EN_binomial_co_1=ENanalysis(y1,X_co,n.rep=10,type_regression="binomial")
save(EN_binomial_co_1,file="../../results/treatment_analysis/EN_binomial_co_1.RData")
#-------------------------------------------------------------------------------
# 2) binomial regression (logistic, y binary) with variable "met_count > 0"
#-------------------------------------------------------------------------------
y2=ifelse(clinical_info$met_count>0,TRUE,FALSE)

EN_binomial_2=ENanalysis(y2,X,n.rep=10,type_regression="binomial")
save(EN_binomial_2,file="../../results/treatment_analysis/EN_binomial_2.RData")
#-------------------------------------------------------------------------------
EN_binomial_co_2=ENanalysis(y2,X_co,n.rep=10,type_regression="binomial")
save(EN_binomial_co_2,file="../../results/treatment_analysis/EN_binomial_co_2.RData")
#-------------------------------------------------------------------------------
# 3) gaussian regression (y = number of metastasis) with variable "met_count"
#-------------------------------------------------------------------------------
y3=clinical_info$met_count

EN_gaussian=ENanalysis(y3,X,n.rep=10,type_regression="gaussian")
save(EN_gaussian,file="../../results/treatment_analysis/EN_gaussian.RData")
#-------------------------------------------------------------------------------
EN_gaussian_co=ENanalysis(y3,X_co,n.rep=10,type_regression="gaussian")
save(EN_gaussian_co,file="../../results/treatment_analysis/EN_gaussian_co.RData")


# RUN PREDICTION WITH DRUGS

# Add drugs as predictors
drug_data <- subset(drug_data, Sample.ID %in% sample_to_keep)
all_drugs <- sort(unique(drug_data$drug))
pred_drugs <- Matrix::t(Matrix::sparseMatrix( i = match(drug_data$drug, all_drugs), 
                                              j = match(drug_data$Sample.ID, sample_to_keep),
                                              x = 1, 
                                              dims = c(length(all_drugs), length(sample_to_keep)),
                                              dimnames = list(all_drugs, sample_to_keep)
))
X_drug <- cbind(X, pred_drugs)
X_co_drug <- cbind(X_co, pred_drugs)

#-------------------------------------------------------------------------------
# 1) binomial regression (logistic, y binary) with variable "is_evidence_of_mets"
#-------------------------------------------------------------------------------
y1=clinical_info$is_evidence_of_mets

EN_binomial_1_drug=ENanalysis(y1,X_drug,n.rep=10,type_regression="binomial")
save(EN_binomial_1_drug,file="../../results/treatment_analysis/EN_binomial_1_drug.RData")
#-------------------------------------------------------------------------------
EN_binomial_co_1_drug=ENanalysis(y1,X_co_drug,n.rep=10,type_regression="binomial")
save(EN_binomial_co_1_drug,file="../../results/treatment_analysis/EN_binomial_co_1_drug.RData")
#-------------------------------------------------------------------------------
# 2) binomial regression (logistic, y binary) with variable "met_count > 0"
#-------------------------------------------------------------------------------
y2=ifelse(clinical_info$met_count>0,TRUE,FALSE)

EN_binomial_2_drug=ENanalysis(y2,X_drug,n.rep=10,type_regression="binomial")
save(EN_binomial_2_drug,file="../../results/treatment_analysis/EN_binomial_2_drug.RData")
#-------------------------------------------------------------------------------
EN_binomial_co_2_drug=ENanalysis(y2,X_co_drug,n.rep=10,type_regression="binomial")
save(EN_binomial_co_2_drug,file="../../results/treatment_analysis/EN_binomial_co_2_drug.RData")
#-------------------------------------------------------------------------------
# 3) gaussian regression (y = number of metastasis) with variable "met_count"
#-------------------------------------------------------------------------------
y3=clinical_info$met_count

EN_gaussian_drug=ENanalysis(y3,X_drug,n.rep=10,type_regression="gaussian")
save(EN_gaussian_drug,file="../../results/treatment_analysis/EN_gaussian_drug.RData")
#-------------------------------------------------------------------------------
EN_gaussian_co_drug=ENanalysis(y3,X_co_drug,n.rep=10,type_regression="gaussian")
save(EN_gaussian_co_drug,file="../../results/treatment_analysis/EN_gaussian_co_drug.RData")


# SAVE SUPPLEMENTARY TABLE 9
openxlsx::write.xlsx(list("predictors_matrix" = as.matrix(X_co_drug), 
                          "predicted_variables" = data.frame( `met_count>0_EN1` = y2,
                                                              met_count_EN2 = y3, row.names = clinical_info$sample_id)
                          ), 
                     file = "../../results/treatment_analysis/supplementary_table_9.xlsx", rowNames=T)











