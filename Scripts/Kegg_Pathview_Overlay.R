############ Libraries ###############
library(pathview)



########################### Set Paths and Create Folders ############################################

# Set base path
base_out <- "/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/Pathway_Overlay"
pathway_id <- "hsa04210"
pathway_name <- "Apoptosis"

# Input files
cd4_file <- "/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/DGE/CD4_BGal_Pos_vs_Neg_DGE.csv"
cd8_file <- "/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/DGE/CD8_BGal_Pos_vs_Neg_DGE.csv"

# Read data
cd4_dge <- read.csv(cd4_file)
cd8_dge <- read.csv(cd8_file)

# Create output dirs
cd4_out <- file.path(base_out, "CD4", pathway_name)
cd8_out <- file.path(base_out, "CD8", pathway_name)
dir.create(cd4_out, recursive = TRUE, showWarnings = FALSE)
dir.create(cd8_out, recursive = TRUE, showWarnings = FALSE)

# Prepare data vectors
cd4_gene_fc <- cd4_dge$log2FoldChange
names(cd4_gene_fc) <- cd4_dge$hgnc_symbol

cd8_gene_fc <- cd8_dge$log2FoldChange
names(cd8_gene_fc) <- cd8_dge$hgnc_symbol


#################################### Run Pathview ###################################################

# Run pathview for CD4
setwd(cd4_out)
pathview(
  gene.data = cd4_gene_fc,
  pathway.id = pathway_id,
  species = "hsa",
  gene.idtype = "SYMBOL",
  out.suffix = "CD4_BGal_Pos_vs_Neg"
)

# Run pathview for CD8
setwd(cd8_out)
pathview(
  gene.data = cd8_gene_fc,
  pathway.id = pathway_id,
  species = "hsa",
  gene.idtype = "SYMBOL",
  out.suffix = "CD8_BGal_Pos_vs_Neg"
)
