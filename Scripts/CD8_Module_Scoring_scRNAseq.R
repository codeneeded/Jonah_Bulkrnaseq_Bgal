## Load libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(clustree)
library(clustifyr)
library(monocle)
library(clusterProfiler)
library(scCustomize)
library(readxl)
library(harmony)
library(DESeq2)
library(vegan)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
options(Seurat.object.assay.version = "v3")



### Conventional CD8+ T cells ###

norm_counts <- readRDS('/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/scRNAseq/conventional_cd8_t_cells/conventional_cd8_rna.rds')
meta <- read.csv('/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/scRNAseq/conventional_cd8_t_cells/conventional_cd8_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/scRNAseq/conventional_cd8_t_cells/conventional_cd8_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.65, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:30, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.8)

DimPlot(data,group.by = 'Cluster_names',repel=T,label=T,label.box = T)
### Create Module Score ############
data <- AddModuleScore(data, features=list(c("CCL3", "SLC4A10", "CDK6", "SOCS2", "CD83", "IL18RAP", "NAMPTP1", "NAMPT", "MIR155HG", "NDFIP2")), name="MS_Tsen", assay='RNA')
data <- AddModuleScore(data, features=list(c("NOG", "SLC16A10", "CACHD1", "CNKSR2","SULT1B1", "LEF1-AS1", "SH3BGRL2", "SH3RF3", "LINC00402", "IL6R")), name="MS_Tctrl", assay='RNA')



### Want to look at it at the sample level

plot_module_by_cluster_sample <- function(seu, modules,
                                          cluster_col = "Cluster_names",
                                          sample_col  = NULL,
                                          point_alpha = 0.6) {
  # infer a sample column if not provided
  if (is.null(sample_col)) {
    candidates <- c("sample","Sample","SampleID","orig.ident","ID","Specimen","Library")
    sample_col <- candidates[candidates %in% colnames(seu@meta.data)][1]
    if (is.na(sample_col)) stop("Please provide sample_col (no common sample column found).")
  }
  
  # pull cluster, sample, and module scores
  vars <- c(cluster_col, sample_col, modules)/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/scRNAseq/
  df <- FetchData(seu, vars = vars)
  
  # long format, average per (cluster, sample, module)
  df_long <- df %>%
    pivot_longer(all_of(modules), names_to = "module", values_to = "score") %>%
    group_by(.data[[cluster_col]], .data[[sample_col]], module) %>%
    summarize(score = mean(score, na.rm = TRUE), .groups = "drop")
  
  # keep cluster order from the object if it's a factor
  if (is.factor(seu@meta.data[[cluster_col]])) {
    df_long[[cluster_col]] <- factor(df_long[[cluster_col]],
                                     levels = levels(seu@meta.data[[cluster_col]]))
  }
  
  # plot: clusters on Y, scores on X
  p <- ggplot(df_long, aes(x = score, y = .data[[cluster_col]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0.15, width = 0, alpha = point_alpha, size = 1) +
    facet_wrap(~ module, scales = "free_x", ncol = 1) +
    labs(x = "Module score (sample mean per cluster)", y = cluster_col) +
    theme_bw()
  
  return(p)
}

setwd('/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/scRNAseq/Module_Score')
# After you ran AddModuleScore:
# data <- AddModuleScore(... name="MS_Tsen", ...)
# data <- AddModuleScore(... name="MS_Tctrl", ...)
# These create columns like MS_Tsen1, MS_Tctrl1

p_A <- plot_module_by_cluster_sample(
  seu        = subset(data, Age_group=='A'),
  modules    = c("MS_Tsen1","MS_Tctrl1"),
  cluster_col= "Cluster_names",
  sample_col = "Donor_id"   # change if your sample column is different
)
ggsave("ModuleScores_byCluster_sampleMean_A.png", p_A, width = 8, height = 10, dpi = 300)

p_E <- plot_module_by_cluster_sample(
  seu        = subset(data, Age_group=='E'),
  modules    = c("MS_Tsen1","MS_Tctrl1"),
  cluster_col= "Cluster_names",
  sample_col = "Donor_id"   # change if your sample column is different
)
ggsave("ModuleScores_byCluster_sampleMean_E.png", p_E, width = 8, height = 10, dpi = 300)

setwd('/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/R_dat')
saveRDS(data, "data.rds", compress = FALSE)  # much bigger file, but very fast

