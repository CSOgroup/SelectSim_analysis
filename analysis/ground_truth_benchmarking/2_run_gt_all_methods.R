#################### Title:

# ground_truth_analysis/2_run_gt_all_methods.R
# Project: "SelectSim_Reviews"

#################### Author:

# Miljan Petrovic, postdoctoral researcher,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# April, 2025.

#################### Description:

# R script to design a ground-truth-including template matrix and run SelectSim
rm(list=ls())
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews")

#################### 

# RUN ALL METHODS ON GIOVANNI'S GROUND TRUTH MODEL

#################### 

# Giovanni's model
parameters <- expand.grid(c(0.95), c(30))
trans_types <- paste0("modelB_", apply(parameters, 1, function(x) paste0("eps",x[1],"_numED",x[2]) ) )

####################
# RUN SelectSim
####################

for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) { 
    
    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))

    lambda_tau <- expand.grid(c(0.3), c(1))
    for (iter_lt in 1:nrow(lambda_tau)) {
      lambda <- lambda_tau[iter_lt,1]
      tau <- lambda_tau[iter_lt,2]
      results <- SelectSim::selectX(  M = run_data_selectSim$M,
                                      sample.class = run_data_selectSim$sample.class,
                                      alteration.class = run_data_selectSim$alteration.class,
                                      n.cores = 1,
                                      min.freq = 5,
                                      n.permut = 1000,
                                      lambda = lambda,
                                      tao = tau,
                                      save.object = FALSE,
                                      verbose = FALSE,
                                      estimate_pairwise = FALSE,
                                      maxFDR = 0.25)$result
      results$name <- apply( results, 1, function(x) paste0(sort(c(x[1], x[2] )), collapse=" - ") )
      results$significance_with_type <- results$type
      results$significance_with_type[results$nFDR2>0.25] <- "none"
      results$observed <- results$overlap
      results$expected <- results$r_overlap
      
      saveRDS(results,
           file=paste0("results/ground_truth_analysis/results_SelectSim_lam",lambda,
                       "_tau",tau,"_ground_truth_",
                       trans_type,gt_type,".rds") )
    }
    print(paste0("Completed ", trans_type, gt_type))
  }
}


####################
# RUN cooccur
####################

for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {
    
    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))
    
    res_cooccur <- cooccur::cooccur(run_data_selectSim$M$M$missense)$results
    res_cooccur$p_lt_gt <- pmax(res_cooccur$p_lt, res_cooccur$p_gt)
    res_cooccur$type <- ifelse(res_cooccur$p_lt_gt > res_cooccur$p_gt, "CO", "ME")
    res_cooccur$p_lt_gt <- res_cooccur$p_lt_gt * ifelse(res_cooccur$type=="CO",1,-1)
    res_cooccur$SFE_1 <- paste0("G", res_cooccur$sp1)
    res_cooccur$SFE_2 <- paste0("G", res_cooccur$sp2)
    res_cooccur$name <- apply( res_cooccur, 1, function(x) paste0(sort(c(x[12], x[13] )), collapse=" - ") )
    res_cooccur$significance_with_type <- "none"
    res_cooccur$significance_with_type[abs(res_cooccur$p_lt_gt)>0.95] <- res_cooccur$type[abs(res_cooccur$p_lt_gt)>0.95]
    res_cooccur$observed <- res_cooccur$obs_cooccur
    res_cooccur$expected <- res_cooccur$exp_cooccur
    
    saveRDS(res_cooccur,
         file=paste0("results/ground_truth_analysis/results_cooccur_ground_truth_",
                     trans_type,gt_type,".rds") )
    print(paste0("Completed ", trans_type, gt_type))
  }
}


####################
# RUN SELECT
####################

for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {
    
    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))
    
    library(select)
    gam <- run_data_selectSim$M$M$missense; gam <- gam[rowSums(gam)>0,colSums(gam)>0]
    res_select <- select::select(M=t(gam), 
                                 sample.class=run_data_selectSim$sample.class,
                                 alteration.class=run_data_selectSim$alteration.class, n.cores = 1, calculate_APC_threshold=F,
                                 save.intermediate.files = F,
                                 verbose = F)
    res_select$significance_with_type <- "none"
    res_select$significance_with_type[res_select$FDR==TRUE] <- res_select$direction[res_select$FDR==TRUE]
    res_select$name <- apply( res_select, 1, function(x) paste0(sort(c(x[1], x[2] )), collapse=" - ") )
    res_select$observed <- res_select$overlap
    res_select$expected <- res_select$r_overlap
    
    saveRDS(res_select,
         file=paste0("results/ground_truth_analysis/results_SELECT_ground_truth_",
                     trans_type,gt_type,".rds") )
    print(paste0("Completed ", trans_type, gt_type))
  }
}


####################
# RUN DISCOVER
####################

for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {
    
    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))
    
    library(discover)
    events <- discover.matrix(run_data_selectSim$M$M$missense)
    result.mutex_ME <- as.data.frame(pairwise.discover.test(events, alternative = "less"))
    result.mutex_CO <- as.data.frame(pairwise.discover.test(events, alternative = "greater"))
    res_discover <- rbind(result.mutex_ME, result.mutex_CO)
    if (nrow(res_discover)==0) {
      res_discover[1,] <- NA
      res_discover$type <- NA
      res_discover$name <- NA
      res_discover$significance_with_type <- "none"
    } else {
      res_discover$type <- c( rep("ME", nrow(result.mutex_ME)), rep("CO", nrow(result.mutex_CO)) )
      res_discover$name <- apply( res_discover, 1, function(x) paste0(sort(c(x[1], x[2] )), collapse=" - ") )
      res_discover$significance_with_type <- "none"
      res_discover$significance_with_type[res_discover$q.value<=0.05] <- res_discover$type[res_discover$q.value<=0.05]
    }
    res_discover$observed <- NA
    res_discover$expected <- NA
    
    saveRDS(res_discover,
            file=paste0("results/ground_truth_analysis/results_DISCOVER_ground_truth_",
                        trans_type,gt_type,".rds") )
    print(paste0("Completed ", trans_type, gt_type))
  }
}


####################
# RUN Coselens
####################
# (requires maf-level data with synonymous and nonsynomnymous mutations so not applicable to our ground-truth gene-level data)
    

####################
# RUN GAMToC (load and results, already ran in Matlab)
####################

# Create config.txt files for GAMToC
for (iter in 1:nrow(parameters)) {
  for (gt_type in c("","_control")) { 
    newdir <- paste0("results/ground_truth_analysis/GAMToC_run_folder/",trans_types[iter],gt_type,"/")
    if (!dir.exists(newdir)) {
      dir.create(newdir, recursive = TRUE)
    }
    conf_text <- write.table(c("save_name = results",
    paste0("patient_mutations = ../../../../data/preprocessed/ground_truth_analysis/gt_maf_",trans_types[iter],gt_type,".txt"),
    "gene_annotation_mat = ../../../../software/gamtoc_code/hg18_exome_compatible.mat",
    "mutation_frequency = 2",
    "algorithm = GREEDY"), 
    file=paste(newdir,"/config_GAMToC_",trans_types[iter],gt_type,".txt",sep=""),
    sep="\n", row.names=F, col.names=F, quote=F)
  }
}

# Load and format results after running in Matlab
for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {
  # for (gt_type in c("")) {
    
    newdir <- paste0("results/ground_truth_analysis/GAMToC_run_folder/",trans_types[iter_trans],gt_type,"/")
    filenames <- list.files(newdir)
    resfile <- filenames[grepl("cytosc",filenames)]
    
    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))
    
    res_gamtoc <- read.table(paste0(newdir,resfile), sep="\t", header=T) 
    res_gamtoc$type <- ifelse(res_gamtoc$corr>0, "CO", "ME")
    res_gamtoc$significance_with_type <- res_gamtoc$type
    res_gamtoc$name <- apply(res_gamtoc, 1, function(x) paste0(c(x[1], x[2] ), collapse=" - ") )
    res_gamtoc$observed <- NA
    res_gamtoc$expected <- NA
    
    saveRDS(res_gamtoc,
            file=paste0("results/ground_truth_analysis/results_GAMToC_ground_truth_",
                        trans_type,gt_type,".rds") )
    print(paste0("Completed ", trans_type, gt_type))
  }
}



#############################
#############################

# RUN WeSME/CO (preprocess GAMs; run in python)

#############################
#############################

# Build input data from GAMs
toy <- read.table("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews/software/wsampling/data/BRCA/BRCA_smut_list.txt")
for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {

    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))

    gam <- run_data_selectSim$M$M$missense

    newdir <- paste0("software/wsampling/data/",trans_types[iter_trans],gt_type,"/")
    if (!dir.exists(newdir)) {
      dir.create(newdir, recursive = TRUE)
    }
    newdir2 <- paste0("software/wsampling/preproc/wrs/",trans_types[iter_trans],gt_type,"/smut/")
    if (!dir.exists(newdir2)) {
      dir.create(newdir2, recursive = TRUE)
    }
    newfilename <- paste0(newdir, trans_types[iter_trans], gt_type, "_smut_list.txt" )
    
    edgelist <- data.frame( gene = rownames(gam),
      sample = apply(gam, 1, function(x) {
      return(paste0( which(x==1), collapse = ","))
    })
    )
    edgelist <- rbind( data.frame(gene = "samples", sample=paste0(colnames(gam),collapse=",")),
                       edgelist)
    data.table::fwrite(edgelist, 
                file=newfilename,
                sep="\t", row.names=F, col.names=F, quote=F)
    
  }
}

# Load and format results after running in Python
for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {
    
    MEs <- read.table(paste0("software/wsampling/results/",trans_types[iter_trans],gt_type,"_smut_me_pvs_0.1.txt"), sep="\t", header=T )
    MEs <- subset(MEs[,c("gene1","gene2","pv..ws.")], pv..ws.<=0.1)
    COs <- read.table(paste0("software/wsampling/results/",trans_types[iter_trans],gt_type,"_smut_co_pvs_0.1.txt"), sep="\t", header=T )
    COs <- subset(COs[,c("gene1","gene2","pv..ws.")], pv..ws.<=0.1)
    MEs$type <- "ME"; COs$type <- "CO"
    MEs$significance_with_type <- "ME"; COs$significance_with_type <- "CO"
    COs$name <- apply(COs, 1, function(x) paste0(c(x[1], x[2] ), collapse=" - ") )
    MEs$name <- apply(MEs, 1, function(x) paste0(c(x[1], x[2] ), collapse=" - ") )
    res_wesme <- rbind(MEs, COs)
    res_wesme$observed <- NA; res_wesme$expected <- NA
    
    saveRDS(res_wesme,
            file=paste0("results/ground_truth_analysis/results_WeSME_ground_truth_",
                        trans_type,gt_type,".rds") )
    print(paste0("Completed ", trans_type, gt_type))
  }
}


####################
# RUN trivial Fisher test
####################

for (iter_trans in 1:length(trans_types)) {  
  trans_type <- trans_types[iter_trans]
  for (gt_type in c("","_control")) {
    
    load(paste0("results/ground_truth_analysis/run_data_gt_",trans_type,gt_type,".RData"))
    gam <- run_data_selectSim$M$M$missense
    
    res_fisher_list <- lapply(as.list(1:nrow(gam)), function(x) {
      tmp_list <- lapply(as.list(1:nrow(gam)), function(y) {
        contig <- table(gam[x,],gam[y,])
        ft_g <- fisher.test(contig, alternative="greater")
        ft_l <- fisher.test(contig, alternative="less")
        df_tmp <- data.frame(gene1 = rownames(gam)[x],
                             gene2 = rownames(gam)[y],
                             pvalue=ifelse(ft_l$p.value<ft_g$p.value, ft_l$p.value, ft_g$p.value),
                             type=ifelse(ft_l$p.value<ft_g$p.value, "ME", "CO")
                           )
      })
      tmp_df <- do.call(rbind, tmp_list)
      return(tmp_df)
    })
    res_fisher <- do.call(rbind, res_fisher_list)
    res_fisher$name <- apply( res_fisher, 1, function(x) paste0(sort(c(x[1], x[2] )), collapse=" - ") )
    res_fisher <- res_fisher[!duplicated(res_fisher$name),]
    res_fisher <- subset(res_fisher, gene1!=gene2)
    res_fisher$qvalue <- p.adjust(res_fisher$pvalue, method="fdr")
    res_fisher$significance_with_type <- res_fisher$type
    res_fisher$significance_with_type[res_fisher$qvalue>0.05] <- "none"
    res_fisher$observed <- NA
    res_fisher$expected <- NA
    
    saveRDS(res_fisher,
            file=paste0("results/ground_truth_analysis/results_Fisher_ground_truth_",
                        trans_type,gt_type,".rds") )
    
    print(paste0("Completed ", trans_type, gt_type))
  }
}





#############################
#############################

##### MERGE ALL RESULTS #####

#############################
#############################

lambda_tau <- expand.grid(c(0.3), c(1))
for (iter_lt in 1:nrow(lambda_tau)) {
  lambda <- lambda_tau[iter_lt,1]
  tau <- lambda_tau[iter_lt,2]
  for (iter_trans in 1:length(trans_types)) {  
    trans_type <- trans_types[iter_trans]
    all_results <- list()
    for (iter_gt_type in 1:2) {
      gt_type <- c("", "_control")[iter_gt_type]
      
      # Load all methods results
      full_results <- list()
      full_results[["SelectSim"]] <- readRDS(paste0("results/ground_truth_analysis/results_SelectSim_lam",lambda,
                                                    "_tau",tau,"_ground_truth_",
                                                    trans_type,gt_type,".rds"))
      full_results[["cooccur"]] <- readRDS(paste0("results/ground_truth_analysis/results_cooccur_ground_truth_",
                                                  trans_type,gt_type,".rds"))
      full_results[["SELECT"]] <- readRDS(paste0("results/ground_truth_analysis/results_SELECT_ground_truth_",
                                                 trans_type,gt_type,".rds"))
      full_results[["DISCOVER"]] <- readRDS(paste0("results/ground_truth_analysis/results_DISCOVER_ground_truth_",
                                                   trans_type,gt_type,".rds"))
      full_results[["GAMToC"]] <- readRDS(paste0("results/ground_truth_analysis/results_GAMToC_ground_truth_",
                                                 trans_type,gt_type,".rds"))
      full_results[["WeSME"]] <- readRDS(paste0("results/ground_truth_analysis/results_WeSME_ground_truth_",
                                                 trans_type,gt_type,".rds"))
      full_results[["Fisher"]] <- readRDS(paste0("results/ground_truth_analysis/results_Fisher_ground_truth_",
                                                trans_type,gt_type,".rds"))
      
      # Merge common columns
      common_columns <- c("name","significance_with_type","observed","expected")
      all_results[[iter_gt_type]] <- rbind( full_results$SelectSim[,common_columns],
                                       full_results$cooccur[,common_columns],
                                       full_results$SELECT[,common_columns],
                                       full_results$DISCOVER[,common_columns],
                                       full_results$GAMToC[,common_columns],
                                       full_results$WeSME[,common_columns],
                                       full_results$Fisher[,common_columns]
      )
      all_results[[iter_gt_type]]$method <- c(rep("SelectSim", nrow(full_results$SelectSim)),
                                         rep("cooccur", nrow(full_results$cooccur)),
                                         rep("SELECT", nrow(full_results$SELECT)),
                                         rep("DISCOVER", nrow(full_results$DISCOVER)),
                                         rep("GAMToC", nrow(full_results$GAMToC)),
                                         rep("WeSME/CO", nrow(full_results$WeSME)),
                                         rep("Fisher", nrow(full_results$Fisher))
      )
      all_results[[iter_gt_type]]$score <- c( full_results$SelectSim$nES,
                                         full_results$cooccur$p_lt_gt,
                                         full_results$SELECT$MI_diff * ifelse(full_results$SELECT$direction=="CO",1,-1),
                                         (full_results$DISCOVER$q.value) * ifelse(full_results$DISCOVER$type=="CO",1,-1),
                                         full_results$GAMToC$corr * ifelse(full_results$GAMToC$type=="CO",1,-1),
                                         (full_results$WeSME$`pv..ws.`) * ifelse(full_results$WeSME$type=="CO",1,-1),
                                         (full_results$Fisher$qvalue) * ifelse(full_results$Fisher$type=="CO",1,-1)
      )
      all_results[[iter_gt_type]]$score_name <- rep( c( "nES",
                                                   "p_lt / -p_gt",
                                                   "Mutual information difference\n(signed by type)",
                                                   "q value (signed by type)",
                                                   "Total correlation",
                                                   "p value (signed by type)",
                                                   "q value (signed by type)"),
                                                   times=c(nrow(full_results$SelectSim),nrow(full_results$cooccur),nrow(full_results$SELECT),nrow(full_results$DISCOVER),nrow(full_results$GAMToC),nrow(full_results$WeSME),nrow(full_results$Fisher)) )
      all_results[[iter_gt_type]]$scores_name <- rep( c( "overlap/r_overlap",
                                                    "obs_cooccur/exp_cooccur",
                                                    "overlap/r_overlap",
                                                    "none",
                                                    "none",
                                                    "none", 
                                                    "none"),
                                                    times=c(nrow(full_results$SelectSim),nrow(full_results$cooccur),nrow(full_results$SELECT),nrow(full_results$DISCOVER),nrow(full_results$GAMToC),nrow(full_results$WeSME),nrow(full_results$Fisher)) )
      all_results[[iter_gt_type]]$data_type <- ifelse( gt_type=="", "GT", "GT_control" )
      
    }
    all_results <- rbind( all_results[[1]], all_results[[2]] )
    
    save(all_results, 
         file = paste0("results/ground_truth_analysis/results_all_methods_ground_truth_",trans_type,
                "_SelSim_lam",lambda,"_tau",tau,
                ".RData"))
    print(paste0("Completed ", trans_type))
  }
}


