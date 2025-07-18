#################### Title:

# normal_mutations/2_run_Lawson.R
# Project: "SelectSim_Reviews"

#################### Author:

# Miljan Petrovic, postdoctoral researcher,
# Ciriello Group, Department of Computational Biology,
# University of Lausanne

# April, 2025.

#################### Description:

# R script to run selectSim on bladder mutation data from Lawson dataset
rm(list=ls())
# setwd("/Users/mpetrov2/mnt/ed2/miljan/SelectSim_Reviews")
library(Matrix)
# source("scripts/rerun_selectSim_new/routines_coincid_test.R")

save_ggplot_custom <- function(plot_name, plot_final, width=40, height=24) {
  ggsave(paste0("", plot_name, ".png"), plot=plot_final, device="png", width=width, height=height, units = "cm")
  return(NA)
}
save_ggplot_custom_pdf <- function(plot_name, plot_final, width=40, height=24) {
  ggsave(paste0("", plot_name, ".pdf"), plot=plot_final, device="pdf", width=width, height=height, units = "cm")
  return(NA)
}
get_ordered_columns <- function(mat, groups) {
  ordered_cols <- c()
  for (grp in levels(groups)) {
    grp_cols <- which(groups == grp)
    if (length(grp_cols) > 1) {
      # Cluster within the group
      dend <- as.dendrogram(hclust(dist(t(mat[, grp_cols]))))
      ordered_grp <- grp_cols[order.dendrogram(dend)]
    } else {
      ordered_grp <- grp_cols
    }
    ordered_cols <- c(ordered_cols, ordered_grp)
  }
  return(ordered_cols)
}
params <- list()
params$cols = c("#66c2a5","#fc8d62","#8da0cb","#d9d57d","#4c6e1b","#7030a0","#946b2d","#0392cf","#f7cac9","#e78ac3","#a6d854","#961203","#AFE1AF","#DC123C","#00468b","#CF9FFF","#BFC0D7")
params$cols4 = c("#AFE1AF","#CF9FFF","#4c6e1b","#7030a0")
params$fsz = 18
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(grid)

# mutset <- "allmut"
mutset <- "cancermut"

#################### 

# Load run data
# load("results/normal_mutations/run_data_Fowler2021.RData")

# Calculate co-incidence & test effect sizes 
# maxFDR <- 0.25
# metric_name <- "overlap"
# r.wobs <- lapply( run_data$nulls, function(x) comut_incidence(x, weighted=TRUE, weights_in=run_data$sample_weights ))
# unwobs <- comut_incidence(run_data$GAM_mut_bin, weighted=FALSE, weights_in=run_data$sample_weights )
# wobs <- comut_incidence(run_data$GAM_mut_bin, weighted=TRUE, weights_in=run_data$sample_weights )
# 
# results <- norm_test_full(unwobs, r.unwobs=NA, wobs, r.wobs, maxFDR, metric_name, run_data$GAM_mut_bin)
# exp.r.nwES <- results$exp.r.nwES
# results <- results$results
# 
# table(results$significance_with_type)

# Run using package SelectSim
runs <- c("normal_samples_all_donors", "normal_samples_normal_donors", "normal_samples_cancer_donors","cancer_samples_cancer_donors")
# for (iter_run in 1:length(runs)) {
#   run <- runs[iter_run]
#   load(paste0("results/normal_mutations/run_data_Lawson2020_",run,".RData"))
#   runtype <- c("T","T","T","C")[iter_run]
#   
#   result_obj <- SelectSim::selectX(  M = run_data_selectSim$M,
#                           sample.class = run_data_selectSim$sample.class,
#                           alteration.class = run_data_selectSim$alteration.class,
#                           n.cores = 1,
#                           min.freq = ceiling(0.01*ncol(run_data_selectSim$M$M$missense)), # 1% of samples
#                           n.permut = 1000,
#                           lambda = 0.3,
#                           tao = 1,
#                           save.object = FALSE,
#                           verbose = FALSE,
#                           estimate_pairwise = FALSE,
#                           maxFDR = 0.25)
#   
#   saveRDS(result_obj, file=paste0("results/normal_mutations/results_Lawson2020_",run,".rds"))
#   
# }

for (iter_run in 1:length(runs)) {
  run <- runs[iter_run]
  load(paste0("results/normal_mutations/run_data_Lawson2020_",run,".RData"))
  runtype <- c("T","T","T","C")[iter_run]
  result_obj <- readRDS(paste0("results/normal_mutations/results_Lawson2020_",run,".rds"))
  
  results_pack <- result_obj$result
  results_pack <- results_pack[order(results_pack$name),]
  results_pack <- results_pack[order(-results_pack$nES),]
  obj_pack <- result_obj$obj

  results_pack$significance_with_type <- results_pack$type
  results_pack$significance_with_type[results_pack$nFDR2>0.25] <- "none"
  plot_ovr_pack <- ggplot( results_pack,
                        aes(x=w_r_overlap, y=w_overlap, 
                            color = significance_with_type,
                        ) )  +
    geom_point(size=1) + 
    scale_color_manual(values=c("ME"=params$cols4[4],"CO"=params$cols4[3],"none"="darkgrey")) + 
    geom_text_repel(aes(label=name), size=5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),) +
    xlab("Mean Simulated Co-incidence") + ylab("Observed Co-incidence") + 
    scale_x_continuous(trans="log1p", breaks = round( exp(seq(log1p(min(results_pack$w_r_overlap)), log1p(max(results_pack$w_r_overlap)), by = 0.5))-1, 1) ) + 
    scale_y_continuous(trans="log1p", breaks = round( exp(seq(log1p(min(results_pack$w_overlap)), log1p(max(results_pack$w_overlap)), by = 0.5))-1, 1) ) + 
    geom_abline(slope=1, intercept=0, color="black", linetype="dashed") +
    guides(color = guide_legend(title = "ED significance")) +
    theme(text = element_text(size = params$fsz), legend.position="top") 
  # plot(plot_ovr_pack)
  
                       
  signif_pairs <- subset( results_pack, name %in% subset(results_pack, significance_with_type %in% c("ME","CO"))$name)
  signif_genes <- unique(c(signif_pairs$SFE_1,signif_pairs$SFE_2))
  plot_mut_count <- ggplot( results_pack,
          aes(x=significance_with_type, y=cum_freq, 
              fill = significance_with_type,
          ) )  +
    geom_boxplot() + #geom_jitter(height=0, width=0.2, size=1, alpha=0.5) +
    scale_fill_manual(values=c("ME"=params$cols4[4],"CO"=params$cols4[3],"none"="grey")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),) +
    xlab(" ") + ylab("Cumulative mutation count") +
    guides(fill = guide_legend(title = "ED significance")) +
    theme(text = element_text(size = params$fsz), legend.position="top") 
  # plot(plot_mut_count)
  
  if (length(signif_genes) > 0) {
  
  ht_genes <- setdiff(signif_genes, c("CDKN1A", "FAT1", "CREBBP"))
  GAM <- pmax(run_data_selectSim$M$M$missense,run_data_selectSim$M$M$truncating)
  tmb_m <- setNames( run_data_selectSim$M$tmb$missense[,"mutation"], run_data_selectSim$M$tmb$missense[,"sample"])
  tmb_t <- setNames( run_data_selectSim$M$tmb$truncating[,"mutation"], run_data_selectSim$M$tmb$truncating[,"sample"])
  blocks <- split(names(run_data_selectSim$sample.class), run_data_selectSim$sample.class)
  blocks <- blocks[order(unlist(lapply(blocks,length)), decreasing = T)]
  ht <- as.matrix(GAM[ht_genes,])
  col2rem <- which(Matrix::colSums(ht)==0)
  ht <- ht[,-col2rem]
  ht <- ht[,order(colnames(ht))]
  group_factor <- factor(run_data_selectSim$sample.class[colnames(ht)], levels = names(blocks))
  ordered_col_indices <- get_ordered_columns(ht, group_factor)
  ht <- ht[,ordered_col_indices]
  colorchange <- which(colSums(ht)==1)
  ht[,colorchange] <- ht[,colorchange]*2
  # how many no mutations samples come from each donor:
  tt <- sapply(unique(run_data_selectSim$sample.class), function(x) { mutcount <- colSums(GAM[ht_genes,which(run_data_selectSim$sample.class==x)]);
  return( c( length(which(mutcount>0)), length(mutcount) )) } ) 
  df_fracs <- data.frame(donor=rep(colnames(tt),2), type=rep(c("samples_with_any_mut","samples_without_mut"),each=ncol(tt)),
                         count = c(tt[1,], tt[2,]-tt[1,]), coloring=c(colnames(tt),rep("A",ncol(tt))),
                         lab = c(paste0(tt[1,],"/",tt[2,]), rep(NA,ncol(tt))) )
  plot_frac <- ggplot( df_fracs,
                       aes(x=donor,y=count, fill=coloring
                       ) )  +
    geom_bar(stat="identity") + geom_label(aes(label=lab), size=3, position = position_stack(vjust = 0.5), show.legend = F) +
    scale_fill_manual(values=setNames(c("grey", RColorBrewer::brewer.pal(12, "Set3"),"skyblue","lightyellow","darkgreen"), c("A",colnames(tt)) ) ) +
    guides(fill = guide_legend(nrow = 3)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    xlab(" ") + ylab("#samples") + ggtitle(paste0("#samples with any mutation within ED genes (colored)\n#samples without mutations (grey)")) +
    theme(text = element_text(size = params$fsz), legend.position="none") + theme(plot.title = element_text(size = 16))
  plot(plot_frac)  
  
  # new "biocportal column order"
  # ...
  
  new_row_order <- order(-rowSums(ht))
  ht <- ht[new_row_order,]
  new_col_order <- c()
  for (ii_block in sort(names(blocks))) {
    
    fullorder <- match(blocks[[ii_block]], colnames(ht))
    if (any(is.na(fullorder))) fullorder <- fullorder[!is.na(fullorder)]
    
    for (iterrow in 1:nrow(ht)) {
      addind <- setdiff( fullorder[which(ht[iterrow,fullorder]>0)], new_col_order )
      new_col_order <- c(new_col_order, addind)
    }
    
  }
  
  ht <- ht[,new_col_order]
  plot_ht <- ComplexHeatmap::Heatmap(ht, 
                    name = "Mutation\nstatus",
                    col = circlize::colorRamp2(c(0,1,2), c("white", "darkred", "darkgreen")),
                    row_names_side = "left",
                    column_title = paste0("Samples (n=", ncol(ht), ", + ", length(col2rem), " w/o mutations)"),
                    row_title = "",
                    show_row_names = F,
                    show_column_names = T,
                    cluster_rows = F, show_row_dend = F,
                    cluster_columns = F, show_column_dend = F,
                    row_gap = grid::unit(3, "pt"),
                    row_split = factor(rownames(ht), levels = rownames(ht)),
                    row_title_gp = grid::gpar(fontsize = 12), 
                    row_title_rot = 0,
                    column_names_gp = grid::gpar(fontsize = 1),
                    heatmap_legend_param = list(title="Mutation\n status", at=c(0,1), labels=c("0","1")),
                    top_annotation = ComplexHeatmap::HeatmapAnnotation(
                      donor = run_data_selectSim$sample.class[colnames(ht)],
                      TMB = ComplexHeatmap::anno_barplot(setNames(tmb_t[colnames(ht)]+tmb_m[colnames(ht)], colnames(ht)),
                                                       gp = grid::gpar(col = "black", fill = "darkblue", color="darkblue"), 
                                                       bar_width = 0.8),
                      annotation_name_side = "left",
                      col = list(donor = setNames(c(RColorBrewer::brewer.pal(12, "Set3"),"skyblue","lightyellow","darkgreen")[1:ncol(tt)], colnames(tt) ))
                    ),
                    show_heatmap_legend = F
  )
  # plot_ht2 <- draw(plot_ht)
  ht_grob <- grid.grabExpr(draw(plot_ht))

  GAM_mut <- readRDS("data/preprocessed/normal_mutations/full_gam_TCGA.rds")
  load("data/raw/TCGA/sample.class_final.RData")
  blca_samples <- split(names(samples.class), samples.class)[["BLCA"]]
  
  load(paste0("results/normal_mutations/GAM_Lawson2020_",run,".RData"))
  signif_genes_sub <- union(intersect(signif_genes, rownames(GAM_mut)), c("TP53","FGFR3","PIK3CA"))
  signif_genes_sub <- intersect(signif_genes_sub, rownames(GAM_full_Lawson))
  gfreq_BLCA <- Matrix::rowSums(GAM_mut[signif_genes_sub,blca_samples]) / length(blca_samples)
  gfreq_pancan <- Matrix::rowSums(GAM_mut[signif_genes_sub,]) / ncol(GAM_mut)
  gfreq_bladder <- Matrix::rowSums(GAM_full_Lawson[signif_genes_sub,]) / ncol(GAM)
  gfreq_bladder_per_pat <- lapply(blocks, function(x) Matrix::rowSums(GAM_full_Lawson[signif_genes_sub,x])/ length(x) )
  gfreq_bladder_per_pat <- data.frame( 
    gene=rep(signif_genes_sub, length(gfreq_bladder_per_pat)),
    mutation_freq = unlist(gfreq_bladder_per_pat),
    donor = rep(names(gfreq_bladder_per_pat), each=length(signif_genes_sub))
  )
  # df_pat <- subset(gfreq_bladder_per_pat, mutation_freq > 0)
  # GAM_pat <- GAM; 
  runtypelabs= c("C"="Cancer", "T"="Healthy")
  # df_freq <- data.frame( gene=c(signif_genes_sub,signif_genes_sub,signif_genes),
  #                        mutation_freq = c(gfreq_BLCA, gfreq_pancan, gfreq_bladder),
  #                        cohort = c(rep("Bladder cancer (TCGA-BLCA)", length(gfreq_BLCA)), rep("All cancers (TCGA)", length(gfreq_pancan)), rep(paste0(runtypelabs[runtype]," bladder (Lawson 2020)"), length(gfreq_bladder)))
  # )
  df_freq <- data.frame( gene=c(signif_genes_sub,signif_genes_sub,signif_genes_sub),
                         mutation_freq = c(gfreq_BLCA, gfreq_pancan, gfreq_bladder),
                         cohort = c(rep("Bladder cancer (TCGA-BLCA)", length(gfreq_BLCA)), rep("All cancers (TCGA)", length(gfreq_pancan)), rep(paste0(runtypelabs[runtype]," bladder (Lawson 2020)"), length(gfreq_bladder)))
  )
  mm <- sum(GAM_mut["FGFR3",blca_samples]) / length(blca_samples)
  df_freq <- rbind(df_freq, data.frame(gene="FGFR3", mutation_freq = mm, cohort = "Bladder cancer (TCGA-BLCA)"))
  plot_mut_freq <- ggplot( subset(df_freq, cohort != "All cancers (TCGA)"),
          aes(x=gene, y=mutation_freq, 
              fill = cohort,
          ) )  +
    geom_bar(stat="identity", position=position_dodge2(width=0.25,padding=0.2)) + 
    # scale_fill_manual(values=c("Bladder cancer (TCGA-BLCA)"=params$cols[7],"All cancers (TCGA)"=params$cols4[12], paste0(runtypelabs[runtype]," bladder (Lawson 2020)")=params$cols[9])) +
    scale_fill_manual(values=c(params$cols[7],params$cols4[12], params$cols[9])) +
    guides(fill = guide_legend(nrow = 3)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),) +
    xlab(" ") + ylab("Gene mutation frequency") + 
    theme(text = element_text(size = params$fsz), legend.position="top") 
  # plot(plot_mut_freq)                       
  plot_mut_freq_per_pat <- ggplot( gfreq_bladder_per_pat,
                           aes(x=gene, y=mutation_freq, 
                               fill = donor,
                           ) )  +
    geom_bar(stat="identity", position=position_dodge2(width=0.25,padding=0.2)) + 
    scale_fill_manual(values=circlize::colorRamp2(c(0,1), c("grey", "darkred"))(seq(0,1,length.out=15))) +
    guides(fill = guide_legend(nrow = 5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlab(" ") + ylab("Gene mutation frequency\n(per donor)") + 
    theme(text = element_text(size = params$fsz), legend.position="top") 
  # plot(plot_mut_freq_per_pat)      
  
  
  # Plot wES normal vs tumor
  results_all_merged <- readRDS("data/raw/selectsim_analysis/analysis_results/catalogue/all_merged_results_v15.rds")
  results_blca <- subset(results_all_merged, Cohort=="TCGA" & Tumor_run=="BLCA")
  new_cols <- paste0(colnames(results_blca), "_BLCA")
  new_rows <- match(results_pack$name, results_blca$name)
  results_pack[which(!is.na(new_rows)),new_cols] <- results_blca[new_rows[!is.na(new_rows)],]
  results_pack$wES_BLCA[is.na(results_pack$wES_BLCA)] <- 0
  results_pack$w_overlap_BLCA[is.na(results_pack$w_overlap_BLCA)] <- 0
  results_pack$significance_with_type_BLCA <- results_pack$type_BLCA
  results_pack$significance_with_type_BLCA[results_pack$nFDR2_BLCA>0.25] <- "none"
  results_pack$significance_with_type_BLCA[is.na(results_pack$significance_with_type_BLCA)] <- "untestable"
  plot_ovr_n2c <- ggplot( subset(results_pack, significance_with_type != "none"),
                           aes(x=wES_BLCA, y=wES, 
                               color = significance_with_type_BLCA,
                               shape = significance_with_type
                           ) )  +
    geom_point(size=4) + 
    scale_color_manual(values=c("ME"=params$cols4[4],"CO"=params$cols4[3],"none"="darkgrey", "untestable"="antiquewhite3")) + 
    geom_text_repel(aes(label=name), size=5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.box = "vertical") +
    xlab("wES (BLCA)") + ylab("wES (Bladder)") + 
    # scale_x_continuous(trans="log1p", breaks = round( exp(seq(log1p(min(results_pack$w_overlap_BLCA)), log1p(max(results_pack$w_overlap_BLCA)), by = 0.5))-1, 1) ) + 
    # scale_y_continuous(trans="log1p", breaks = round( exp(seq(log1p(min(results_pack$w_overlap)), log1p(max(results_pack$w_overlap)), by = 0.5))-1, 1) ) + 
    geom_abline(slope=1, intercept=0, color="black", linetype="dashed") +
    guides(color = guide_legend(title = "ED significance in tumor"), shape = guide_legend(title = "ED significance in healthy")) +
    theme(text = element_text(size = params$fsz), legend.position="top") 
  # plot(plot_ovr_n2c)
  
  
  plot_final <- grid.arrange(
    grobs = list(plot_ovr_pack, plot_mut_count, plot_mut_freq, ht_grob, plot_mut_freq_per_pat, plot_ovr_n2c,plot_frac),
    layout_matrix = rbind(c(1, 2, 3, 5, 6),
                          c(1, 4, 4, 7, 7)) )
  } else {
    
    plot_final <- plot_ovr_pack
  }
  save_ggplot_custom(paste0("figures/normal_mutations/Lawson2020_summary_",run), plot_final, width=56, height=24)
  save_ggplot_custom_pdf(paste0("figures/normal_mutations/Lawson2020_summary_",run), plot_final, width=56, height=24)
  
}
  


















# compare for muscle-invasiveness
for (iter_plot in 1:2) {
  
  maf <- readRDS("data/preprocessed/normal_mutations/lawson2020_maf.rds")
  if (iter_plot == 1) maf <- subset(maf, histology=="Urothelium")
  if (iter_plot == 2) maf <- subset(maf, histology!="Urothelium")
  sub2samp <- split(maf$sample, maf$subject)
  sub_sizes <- unlist(lapply(sub2samp, function(x) return(length(unique(x)))))
  df <- expand.grid(sort(unique(maf$gene)), sort(unique(maf$subject)) ); colnames(df) <- c("gene", "subject")
  df$num_samples <- NA
  df$num_muts <- NA
  for (iter in 1:nrow(df)) {
    tmp_maf <- subset(maf, subject == df$subject[iter] & gene == df$gene[iter])
    df$num_samples[iter] <- length(unique(tmp_maf$sample))
    df$altfreq[iter] <- df$num_samples[iter] / sub_sizes[df$subject[iter]]
    df$num_muts[iter] <- length(unique(tmp_maf$AA_change))
    
    # print(iter)
  }
  df$muscle_invasive <- maf$subject_type[ match(df$subject, maf$subject)]
  
  
  # plot heatmap of donors and mutations
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. Prepare matrices from your dataframe (df)
  genes <- sort(unique(df$gene))
  subjects <- sort(unique(df$subject))
  mat_altfreq <- as.matrix(sparseMatrix(i = match(df$gene, genes), 
                              j = match(df$subject, subjects),
                              x = as.numeric(df$altfreq), 
                              dimnames = list(genes, subjects)))
  mat_num_samples <- as.matrix(sparseMatrix(i = match(df$gene, genes), 
                              j = match(df$subject, subjects),
                              x = df$num_samples, 
                              dimnames = list(genes, subjects)))
  mat_num_muts <- as.matrix(sparseMatrix(i = match(df$gene, genes), 
                              j = match(df$subject, subjects),
                              x = df$num_muts, 
                              dimnames = list(genes, subjects)))
  
  
  # 2. Define color mapping for 'altfreq'
  col_fun <- colorRamp2(
    breaks = c(min(mat_altfreq, na.rm = TRUE), max(mat_altfreq, na.rm = TRUE)),
    colors = c(params$cols[15], params$cols[14]) # Customize colors as needed
  )
  if (iter_plot == 1) cairo_pdf(paste0("figures/normal_mutations/Lawson2020_gene_alt_freq_vs_all_donors_normal.pdf"), width = 15, height = 15)
  if (iter_plot == 2) cairo_pdf(paste0("figures/normal_mutations/Lawson2020_gene_alt_freq_vs_all_donors_cancer.pdf"), width = 15, height = 15)
  if (iter_plot == 1) {title_ht <- "fraction of normal samples\nwith a mutation"; fontsize <- 5; fontsize_row <- 5}
  if (iter_plot == 2) {title_ht <- "fraction of tumor samples\nwith a mutation"; fontsize <- 8; fontsize_row <- 10}
  ht_plot <- Heatmap(
    mat_altfreq,
    name = title_ht, # Legend title
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        label = paste0(mat_num_samples[i, j], "/", mat_num_muts[i, j]),
        x = x, y = y,
        gp = gpar(fontsize = fontsize) # Adjust font size as needed
      )
    },
    show_row_dend = F, show_column_dend = F,
    row_title = " ",
    column_title = "<<#samples w/ mut.>> / <<#distinct mut.>>",
    row_names_side = "left",
    row_names_gp = gpar(fontsize = fontsize_row),
    column_names_side = "bottom",
    show_heatmap_legend = TRUE,
    top_annotation = HeatmapAnnotation(
      muscle_invasiveness = setNames( df$muscle_invasive[match(subjects, df$subject)], subjects )
    )
  )
  ht_plot <- draw(ht_plot,padding = unit(c(2, 2, 2, 2), "cm") )
  dev.off()
  
}


