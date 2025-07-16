#################### Title:

# ground_truth_analysis/1_design_ground_truth_modelB.R
# Project: "SelectSim_Reviews"

#################### Author:

# Giovanni Ciriello, Miljan Petrovic,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# April, 2025.

#################### Description:

# R script to design a ground-truth 
rm(list=ls())
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews")

#################### 

set.seed(42)

# initialize parameters

ngenes = 100
nsamples = 1000
# ngenes = 50
# nsamples = 500

# min and max tumor mutation burdens
min_tmb = 5
max_tmb = 30
# min_tmb = 1
# max_tmb = 10

# min and max gene alteration frequencies
min_gfr = 0.05
max_gfr = 0.2
# min_gfr = 0.05
# max_gfr = 0.2

# eps sets the influence of CO/ME relationships on mutation selection
# eps = 0.5 will increase or decrease the probability of acquiring a given mutation by 50%
eps_s <- c(0.75, 0.85, 0.95)
numED_s <- c(20,30,40)
# eps_s <- c(0.75, 0.85, 0.95)
# numED_s <- c(10,20)
# eps_s <- c(0.75)
# numED_s <- c(10)
for ( eps in eps_s) {
  for (numED in numED_s) {
    
    # initialize variables
    
    gam = matrix(0, nrow = ngenes, ncol = nsamples)
    
    tmb = sample(min_tmb:max_tmb, size = nsamples, replace = TRUE)
    # tmb = sample(min_tmb:max_tmb, size = nsamples, replace = TRUE, prob = 1/((min_tmb:max_tmb)^1.05) )
    sum(tmb) / (ngenes*nsamples)
    gfr = sample(seq(from = min_gfr, to = max_gfr, by = 0.01), size = ngenes, replace = TRUE)
    # gfr = sample(seq(from = min_gfr, to = max_gfr, by = 0.01), size = ngenes, replace = TRUE, prob = 1/(seq(from = min_gfr, to = max_gfr, by = 0.01)^3) )
    
    # here I sample what will be our true evolutionary dependencies
    # 20 ED with 10 ME and 10 CO
    ed = matrix(0, nrow = nrow(gam), ncol = nrow(gam))
    rownames(ed) = sapply(1:nrow(ed), function(x) paste('G',x,sep=''))
    colnames(ed) = sapply(1:nrow(ed), function(x) paste('G',x,sep=''))
    
    pairs = matrix(, nrow = numED, ncol = 3)
    
    prob = rep(1/ngenes, nrow(ed))
    for(i in 1:nrow(pairs)){
      
      value = 1
      if(i%%2 == 0)
        value = -1
      
      genes = sample(1:nrow(ed), 2, replace = F, prob = prob)
      ed[genes[1],genes[2]] = value
      ed[genes[2],genes[1]] = value
      
      # by setting the probabilities of the genes that were already sampled to zero, we avoid creating paths or loop which may lead to transitive effects
      prob[genes[1]] = 0
      prob[genes[2]] = 0
      prob[prob != 0] = 1/sum(prob != 0)
      
      pairs[i, ] = c(rownames(ed)[genes[1]],rownames(ed)[genes[2]],value)
      
    }
    
    
    # here I simulate a GAM, one sample at a time
    
    rownames(gam) = rownames(ed)
    colnames(gam) = sapply(1:ncol(gam), function(x) paste('S',x,sep=''))
    
    for(i in 1:ncol(gam)){
      
      # init variables
      prob = gfr
      nmut = 0
      
      while(nmut < tmb[i]){
        
        position = sample(1:nrow(gam), 1, replace=F, prob=prob)
        gam[position,i] = 1
        nmut = nmut + 1
        
        # update prob vector based on ED
        # I also set the probability of mutating twice the same gene equal to zero
        prob = prob + prob*(eps*ed[position,])
        prob[position] = 0

        # this simply renormalize the prob vector to make it sum to 1
        prob[prob != 0] = prob[prob != 0]/sum(prob[prob != 0])
        
      }
      
    }
    
    # create list of control gam without ED
    
    ctrl = list()
    
    for(k in 1:100){
      
      gam_control = matrix(0, nrow = ngenes, ncol = nsamples)
      rownames(gam_control) = rownames(gam)
      colnames(gam_control) = colnames(gam)
      
      for(i in 1:ncol(gam_control)){
        
        # init variables
        prob = gfr
        nmut = 0
        
        while(nmut < tmb[i]){
          
          position = sample(1:nrow(gam_control), 1, replace=F, prob=prob)
          gam_control[position,i] = 1
          nmut = nmut + 1
          
          # update prob vector	
          prob[position] = 0
          prob[prob != 0] = prob[prob != 0]/sum(prob[prob != 0])
          
        }
        
      }
      
      ctrl[[k]] = gam_control
      
    }
    
    genes <- rownames(gam); samples <- colnames(gam)
    
    # Build data for SelectSim (gam)
    run_data_selectSim <- list()
    run_data_selectSim$M <- list()
    run_data_selectSim$M$M <- list()
    run_data_selectSim$M$M$missense <- as.matrix(gam)
    run_data_selectSim$M$M$truncating <- matrix(0, nrow=nrow(gam), ncol=ncol(gam), dimnames=dimnames(gam) )
    run_data_selectSim$M$tmb <- list()
    run_data_selectSim$M$tmb$missense <- data.frame( sample = samples, mutation = colSums(gam), row.names = samples )
    # run_data_selectSim$M$tmb$truncating <- data.frame( sample = samples, mutation = rep(1, length(samples)), row.names = samples )
    run_data_selectSim$M$tmb$truncating <- data.frame( sample = samples, mutation = rep(0.00001, length(samples)), row.names = samples )
    run_data_selectSim$sample.class <- setNames(rep("TUMOR_TYPE",length(samples)), samples)
    run_data_selectSim$alteration.class <- setNames(rep("MUT", length(genes)), genes)
    save(run_data_selectSim, pairs,
         file = paste0("results/ground_truth_analysis/run_data_gt_modelB_eps",eps,"_numED",numED,".RData"))
    
    # Build data for GAMToC
    ii <- which(gam!=0, arr.ind=T)
    write.table( data.frame(case=colnames(gam)[ii[,2]], gene=rownames(gam)[ii[,1]] ),
                 file=paste0("data/preprocessed/ground_truth_analysis/gt_maf_modelB_eps",eps,"_numED",numED,".txt"), sep="\t", row.names=F, col.names=T, quote=F )
    
    # Build data for SelectSim (gam_control)
    run_data_selectSim <- list()
    run_data_selectSim$M <- list()
    run_data_selectSim$M$M <- list()
    run_data_selectSim$M$M$missense <- as.matrix(gam_control)
    run_data_selectSim$M$M$truncating <- matrix(0, nrow=nrow(gam), ncol=ncol(gam), dimnames=dimnames(gam) )
    run_data_selectSim$M$tmb <- list()
    run_data_selectSim$M$tmb$missense <- data.frame( sample = samples, mutation = colSums(gam_control), row.names = samples )
    run_data_selectSim$M$tmb$truncating <- data.frame( sample = samples, mutation = rep(0.00001, length(samples)), row.names = samples )
    run_data_selectSim$sample.class <- setNames(rep("TUMOR_TYPE",length(samples)), samples)
    run_data_selectSim$alteration.class <- setNames(rep("MUT", length(genes)), genes)
    save(run_data_selectSim, ctrl, file = paste0("results/ground_truth_analysis/run_data_gt_modelB_eps",eps,"_numED",numED,"_control.RData"))
    
    # Build data for GAMToC
    ii <- which(gam_control!=0, arr.ind=T)
    write.table( data.frame(case=colnames(gam_control)[ii[,2]], gene=rownames(gam_control)[ii[,1]] ),
                 file=paste0("data/preprocessed/ground_truth_analysis/gt_maf_modelB_eps",eps,"_numED",numED,"_control.txt"), sep="\t", row.names=F, col.names=T, quote=F )
    
    print(c(eps,numED))
  }  
}



