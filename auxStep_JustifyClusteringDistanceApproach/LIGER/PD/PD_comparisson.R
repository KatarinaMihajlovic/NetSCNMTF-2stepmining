library(rliger)
library(Matrix) 
library(patchwork)
library(Seurat)
library(glue)

days <- list("IPSCs", "D06", "D15", "D21")
days_d <- c("IPSCs"="D0", "D06"="D6", "D15"="D15", "D21"="D21")
for (day in days) {
  day
  ctrl_dge <- read.csv(glue("input/E_Control_{day}.csv"), row.names = 1, header= TRUE);
  ctrl_dge <- data.matrix(ctrl_dge)
  ctrl_dge <- as(ctrl_dge, "dgTMatrix")
  stim_dge <- read.csv(glue("input/E_PINK1_{day}.csv"), row.names = 1, header= TRUE);
  stim_dge <- data.matrix(stim_dge)
  stim_dge <- as(stim_dge, "dgTMatrix")
  # matrix_list <- read10X(sample.dirs =c("10x_ctrl_outs", "10x_stim_outs"), sample.names = c("ctrl", "stim"), merge = F);
  ifnb_liger <- createLiger(list(ctrl = ctrl_dge, stim = stim_dge))
  
  # ifnb_liger <- normalize(ifnb_liger)
  ifnb_liger@norm.data <- ifnb_liger@raw.data 
  ifnb_liger <- selectGenes(ifnb_liger)
  ifnb_liger <- scaleNotCenter(ifnb_liger)
  
  ifnb_liger <- optimizeALS(ifnb_liger, k = 100)
  ifnb_liger <- quantile_norm(ifnb_liger)
  
  ifnb_liger <- runUMAP(ifnb_liger, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
  
  all.plots <- plotByDatasetAndCluster(ifnb_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
  all.plots[[1]] + all.plots[[2]]
  
  #gene_loadings <- plotGeneLoadings(ifnb_liger, do.spec.plot = FALSE, return.plots = TRUE)
  #gene_loadings[[1]]
  
  
  #datsaets
  datasets.results <- runWilcoxon(ifnb_liger, compare.method = "datasets")
  head(datasets.results)
  
  datasets.results <- datasets.results[datasets.results$padj < 0.05,]
  
  datasets.results <- datasets.results[datasets.results$logFC > 3,]
  
  wilcoxon.datasets <- datasets.results[order(datasets.results$padj), ]
  markers <- wilcoxon.datasets[1:20, ]
  head(markers)
  
  #ISG15 <- plotGene(ifnb_liger, "ISG15", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
  #ISG15[[1]] + ISG15[[2]]
  write.csv(wilcoxon.datasets,glue("output/SS_PDpreds_{days_d[day]}.csv"), row.names = FALSE)
}




