# Author: Debora Sesisa
# EN MSK data

source("/Users/dsesia/Desktop/Analysis per topic/EN_Arvind/LUAD/code/EN_function.R")

load("/Users/dsesia/Desktop/Analysis per topic/EN_Arvind/all_MSK/data/PM_comutations.RData")
load("/Users/dsesia/Desktop/Analysis per topic/EN_Arvind/all_MSK/data/full_gam.RData")
clinical_info=read.delim("/Users/dsesia/Desktop/Analysis per topic/EN_Arvind/all_MSK/data/msk_full_clinical_table.txt")

# code for server
source("functions/EN_function.R")

load("data/PM_comutations.RData")
load("data/full_gam.RData")
clinical_info=read.delim("data/msk_full_clinical_table.txt")

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

#-------------------------------------------------------------------------------
# 1) binomial regression (logistic, y binary) with variable "is_evidence_of_mets"
#-------------------------------------------------------------------------------
y1=clinical_info$is_evidence_of_mets

EN_binomial_1=ENanalysis(y1,X,n.rep=10,type_regression="binomial")
save(EN_binomial_1,file="results_data/EN_binomial_1.RData")
#-------------------------------------------------------------------------------
EN_binomial_co_1=ENanalysis(y1,X_co,n.rep=10,type_regression="binomial")
save(EN_binomial_co_1,file="results_data/EN_binomial_co_1.RData")
#-------------------------------------------------------------------------------
# 2) binomial regression (logistic, y binary) with variable "met_count > 0"
#-------------------------------------------------------------------------------
y2=ifelse(clinical_info$met_count>0,TRUE,FALSE)

EN_binomial_2=ENanalysis(y2,X,n.rep=10,type_regression="binomial")
save(EN_binomial_2,file="results_data/EN_binomial_2.RData")
#-------------------------------------------------------------------------------
EN_binomial_co_2=ENanalysis(y2,X_co,n.rep=10,type_regression="binomial")
save(EN_binomial_co_2,file="results_data/EN_binomial_co_2.RData")
#-------------------------------------------------------------------------------
# 3) gaussian regression (y = number of metastasis) with variable "met_count"
#-------------------------------------------------------------------------------
y3=clinical_info$met_count

EN_gaussian=ENanalysis(y3,X,n.rep=10,type_regression="gaussian")
save(EN_gaussian,file="results_data/EN_gaussian.RData")
#-------------------------------------------------------------------------------
EN_gaussian_co=ENanalysis(y3,X_co,n.rep=10,type_regression="gaussian")
save(EN_gaussian_co,file="results_data/EN_gaussian_co.RData")
