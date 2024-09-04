library('dplyr')
library('tidyr')
library('ggplot2')
library('ggpubr')
library('tictoc')
library('coselens')
library('foreach')
library('doParallel')

# Note: coselens package used over here is realse v1

# Functions to run coselens
coselns_test <- function(genes,maf,gene_list){
    require('dplyr')
    require('tidyr')
    `%notin%` <- Negate(`%in%`)
    print(genes)
    mut_samples<-unique((maf %>% filter(gene==genes))$sampleID)
    group_1 <- maf %>% filter(sampleID %in% mut_samples) %>% select('sampleID','chr','pos','ref','alt')
    group_2 <- maf %>% filter(sampleID %notin% mut_samples) %>% select('sampleID','chr','pos','ref','alt')
    group_1$pos <- as.numeric(group_1$pos)
    group_2$pos <- as.numeric(group_2$pos)
    result <- coselens(group_1, group_2,gene_list)
    if(nrow(result)>=1){
        temp <-result %>% mutate('DeltaNd'=num.drivers.group1-num.drivers.group2) %>% mutate('type'=case_when(DeltaNd<0 ~ 'ME',DeltaNd>0 ~ 'CO')) %>% 
                    mutate('split_gene'=rep(genes,n())) %>% mutate('pair'=paste(split_gene,gene_name,sep=" - "))
        return(temp)
    }
    else{
            return(result)
    }


}

# Gene to consider
gene_list <-readRDS('/mnt/ptemp/arvind/catalouge_work/gene_used.rds')
# Maf file
luad_maf_data<- readRDS(file='/mnt/ptemp/arvind/tool_comaprision/data/luad_maf_coselns_data.rds')

# Gene to consider
gene_to_consider <- c()
i=1
for(genes in gene_list[[1]]){
        mut_samples<-unique((luad_maf_data %>% filter(gene==genes))$sampleID)
        if(length(mut_samples)>=25){
            gene_to_consider<-c(genes,gene_to_consider)
            i=i+1
        }
}

print(length(gene_to_consider))
print(gene_to_consider[1:5])

registerDoParallel(5)  # use multicore, set to the number of our cores

tic('Total Time:')
pan_result<-foreach(i = c(1:length(gene_to_consider)),.combine = rbind) %dopar% {
  coselns_test(genes = gene_to_consider[[i]],maf = luad_maf_data,gene_list= unlist(gene_list[[1]]))
}
toc()

print(dim(pan_result))
saveRDS(pan_result,file='/mnt/ptemp/arvind/tool_comaprision/results/coselens/coselns_luad_all_results.rds')
#coselns_test(genes = gene_to_consider[[1]],maf = luad_maf_data,gene_list= unlist(gene_list[[1]]))


# Smapling results
sample_list1 <- readRDS('/mnt/ptemp/arvind/tool_comaprision/data/sampling_list_luad.rds')
sample_list2 <- readRDS('/mnt/ptemp/arvind/tool_comaprision/data/sampling_list_luad_2.rds')
sample_list<-c(sample_list1[11:15],sample_list2[11:15])

registerDoParallel(20)  # use multicore, set to the number of our cores
result_vec<-list()
k=10
tic('Total Time:')
for(i in c(1:10)){
    sample_data <- luad_maf_data %>% filter(sampleID %in% names(sample_list[[i]]))
    tic('Run Time:')
        pan_result<-foreach(i = c(1:length(gene_to_consider)),.combine = rbind) %dopar% {
            coselns_test(genes = gene_to_consider[[i]],maf = sample_data,gene_list= unlist(gene_list[[1]]))
        }
    saveRDS(pan_result,file=paste('/mnt/ptemp/arvind/tool_comaprision/results/coselens/luad_sampling_run_',k,'_v2.rds',sep=""))
    toc()
    result_vec[[k]]<-pan_result
    k=k+1
}
tic('Total Time:')
saveRDS(result_vec,file='/mnt/ptemp/arvind/tool_comaprision/results/coselens/luad_sampling_results_v2.rds')