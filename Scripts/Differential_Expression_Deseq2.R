# Libraries 

library(DESeq2)
library(tidyverse)
library(biomaRt)
library(EnhancedVolcano)
library(VennDiagram)
library(pheatmap)
library(dplyr)

# Load count data (raw featureCounts output)
counts_df <- read.csv("/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/Feature Counts/raw_combined_featureCounts.csv", row.names = 1, check.names = FALSE)


########################################################### Data Cleaning ###################################################
# Load metadata
metadata <- read.csv("/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/Bulk_RNA-Seq_Metadata.csv")
colnames(metadata) <- make.names(colnames(metadata))  # Sanitize column names

# Rename 'B.Gal' to 'B_Gal'
colnames(metadata)[colnames(metadata) == "B.Gal"] <- "B_Gal"
# Standardize sample names to lowercase
metadata$File_Name <- tolower(metadata$File_Name)

# Filter metadata and count matrix to matching samples
metadata <- metadata %>% filter(File_Name %in% colnames(counts_df))
counts_df <- counts_df[, metadata$File_Name]

# Set rownames so DESeq2 can align properly
rownames(metadata) <- metadata$File_Name


# After loading and aligning metadata
metadata$ID <- factor(metadata$ID)
metadata$B_Gal <- factor(metadata$B_Gal)

# Clean and standardize B_Gal labels
metadata$B_Gal <- recode(metadata$B_Gal,
                         "Pos" = "Positive",
                         "Neg" = "Negative")
################################################# Differential Gene Expression #############################################################

# Set output directory
dge_output_dir <- "/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/DGE"
if (!dir.exists(dge_output_dir)) dir.create(dge_output_dir, recursive = TRUE)

# Set up biomaRt for human gene annotation
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# ---------- CD4 ----------
metadata_cd4 <- metadata %>% filter(Cell_Type == "CD4")
counts_cd4 <- counts_df[, metadata_cd4$File_Name]
rownames(metadata_cd4) <- metadata_cd4$File_Name

dds_cd4 <- DESeqDataSetFromMatrix(
  countData = counts_cd4,
  colData = metadata_cd4,
  design = ~ ID + B_Gal
)
dds_cd4 <- DESeq(dds_cd4)

res_cd4 <- results(dds_cd4, contrast = c("B_Gal", "Positive", "Negative"))
res_cd4_df <- as.data.frame(res_cd4) %>%
  rownames_to_column("Ensembl_ID") %>%
  arrange(padj)

# Annotate
symbols_cd4 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_cd4_df$Ensembl_ID,
  mart = ensembl
)

# Merge and save
res_cd4_annotated <- res_cd4_df %>%
  left_join(symbols_cd4, by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  relocate(hgnc_symbol, .after = Ensembl_ID)

write.csv(res_cd4_annotated, file = file.path(dge_output_dir, "CD4_BGal_Pos_vs_Neg_DGE.csv"), row.names = FALSE)


# ---------- CD8 ----------
metadata_cd8 <- metadata %>% filter(Cell_Type == "CD8")
counts_cd8 <- counts_df[, metadata_cd8$File_Name]
rownames(metadata_cd8) <- metadata_cd8$File_Name

dds_cd8 <- DESeqDataSetFromMatrix(
  countData = counts_cd8,
  colData = metadata_cd8,
  design = ~ ID + B_Gal
)
dds_cd8 <- DESeq(dds_cd8)

res_cd8 <- results(dds_cd8, contrast = c("B_Gal", "Positive", "Negative"))
res_cd8_df <- as.data.frame(res_cd8) %>%
  rownames_to_column("Ensembl_ID") %>%
  arrange(padj)

# Annotate
symbols_cd8 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_cd8_df$Ensembl_ID,
  mart = ensembl
)

# Merge and save
res_cd8_annotated <- res_cd8_df %>%
  left_join(symbols_cd8, by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  relocate(hgnc_symbol, .after = Ensembl_ID)

write.csv(res_cd8_annotated, file = file.path(dge_output_dir, "CD8_BGal_Pos_vs_Neg_DGE.csv"), row.names = FALSE)


############################################################ PCA ###########################################################


# 1. Add a combined group column for coloring
metadata$Group <- paste0(metadata$Cell_Type, "_", metadata$B_Gal)

# 2. Subset to samples with count data
counts_pca <- counts_df[, metadata$File_Name]
rownames(metadata) <- metadata$File_Name

# 3. Create DESeqDataSet for PCA only (no need to run DESeq)
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_pca,
  colData = metadata,
  design = ~ Group
)

# 4. rlog transform for PCA
rld_all <- rlog(dds_all, blind = TRUE)

# 5. Extract PCA coordinates
pca_data <- plotPCA(rld_all, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# 6. Plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(type = "norm", level = 0.95, size = 1.2, linetype = "dashed", alpha = 0.4) +
  geom_point(size = 4, alpha = 0.85) +
  labs(
    title = "PCA of All Samples by Cell Type and B-Gal Status",
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    color = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set1")

ggsave(
  filename = "/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/DGE/PCA_CellType_BGal.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

################################ Volcano Plot ####################################################

# ---------- CD4 Volcano ----------
# Filter for genes with symbol
cd4_deg <- res_cd4_annotated %>% filter(!is.na(hgnc_symbol))

# Custom coloring
cd4_keyvals <- ifelse(
  abs(cd4_deg$log2FoldChange) > 1.5 & cd4_deg$padj < 0.01, "#CD0BBC",
  ifelse(cd4_deg$padj < 0.01, "#28E2E5", 'gray30')
)
cd4_keyvals[is.na(cd4_keyvals)] <- 'gray30'
names(cd4_keyvals)[cd4_keyvals == 'gray30'] <- 'NS'
names(cd4_keyvals)[cd4_keyvals == '#28E2E5'] <- 'adj(p-value) < 0.01'
names(cd4_keyvals)[cd4_keyvals == '#CD0BBC'] <- 'FC > 1.5'

# Plot
vp_cd4 <- EnhancedVolcano(cd4_deg,
                          lab = cd4_deg$hgnc_symbol,
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          pCutoffCol = 'padj',
                          xlab = bquote(~Log[2]~ 'fold change'),
                          pCutoff = 0.01,
                          FCcutoff = 1.5,
                          pointSize = 4.0,
                          labSize = 3.0,
                          labCol = 'black',
                          labFace = 'bold',
                          colAlpha = 4/5,
                          legendPosition = 'right',
                          legendLabSize = 14,
                          legendIconSize = 4.0,
                          drawConnectors = TRUE,
                          widthConnectors = 1.0,
                          colConnectors = 'black',
                          title = 'CD4: B-Gal+ vs B-Gal-',
                          subtitle = 'Differential Gene Expression',
                          colCustom = cd4_keyvals
)

# Save
ggsave(
  filename = file.path(dge_output_dir, "CD4_BGal_Volcano.png"),
  plot = vp_cd4 + guides(color = guide_legend(reverse = TRUE)),
  dpi = 500,
  width = 10,
  height = 7
)


# ---------- CD8 Volcano ----------
cd8_deg <- res_cd8_annotated %>% filter(!is.na(hgnc_symbol))

cd8_keyvals <- ifelse(
  abs(cd8_deg$log2FoldChange) > 1.5 & cd8_deg$padj < 0.01, "#CD0BBC",
  ifelse(cd8_deg$padj < 0.01, "#28E2E5", 'gray30')
)
cd8_keyvals[is.na(cd8_keyvals)] <- 'gray30'
names(cd8_keyvals)[cd8_keyvals == 'gray30'] <- 'NS'
names(cd8_keyvals)[cd8_keyvals == '#28E2E5'] <- 'adj(p-value) < 0.01'
names(cd8_keyvals)[cd8_keyvals == '#CD0BBC'] <- 'FC > 1.5'

vp_cd8 <- EnhancedVolcano(cd8_deg,
                          lab = cd8_deg$hgnc_symbol,
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          pCutoffCol = 'padj',
                          xlab = bquote(~Log[2]~ 'fold change'),
                          pCutoff = 0.01,
                          FCcutoff = 1.5,
                          pointSize = 4.0,
                          labSize = 3.0,
                          labCol = 'black',
                          labFace = 'bold',
                          colAlpha = 4/5,
                          legendPosition = 'right',
                          legendLabSize = 14,
                          legendIconSize = 4.0,
                          drawConnectors = TRUE,
                          widthConnectors = 1.0,
                          colConnectors = 'black',
                          title = 'CD8: B-Gal+ vs B-Gal-',
                          subtitle = 'Differential Gene Expression',
                          colCustom = cd8_keyvals
)

ggsave(
  filename = file.path(dge_output_dir, "CD8_BGal_Volcano.png"),
  plot = vp_cd8 + guides(color = guide_legend(reverse = TRUE)),
  dpi = 500,
  width = 10,
  height = 7
)

######################################### Venn Diagram - CD4 CD8 Bgal interactions ############################
# ---- 1. Filter and rank CD4 DEGs ----
cd4_top <- res_cd4_annotated %>%
  filter(!is.na(hgnc_symbol), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 200)

# ---- 2. Filter and rank CD8 DEGs ----
cd8_top <- res_cd8_annotated %>%
  filter(!is.na(hgnc_symbol), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 200)

# ---- 3. Create gene symbol sets ----
genes_cd4 <- unique(cd4_top$hgnc_symbol)
genes_cd8 <- unique(cd8_top$hgnc_symbol)

# ---- 4. Intersecting genes ----
shared_genes <- intersect(genes_cd4, genes_cd8)

# ---- 5. Save intersecting genes ----
write.csv(shared_genes, file = file.path(dge_output_dir, "Shared_DEGs_CD4_CD8_BGal_Pos.csv"), row.names = FALSE)

# ---- 6. Venn Diagram (adjusted layout) ----
venn.plot <- venn.diagram(
  x = list(
    "CD4 B-Gal+" = genes_cd4,
    "CD8 B-Gal+" = genes_cd8
  ),
  filename = NULL,
  fill = c("#E41A1C", "#377EB8"),  # red + blue
  alpha = 0.6,
  cex = 2,
  cat.cex = 2,
  cat.pos = c(0, 0),         # Center category labels horizontally
  cat.dist = 0.03,           # Bring category labels closer
  margin = 0.05,             # Reduce white space around plot
  lwd = 2,
  col = "black",
  main = "Overlap of Top DEGs in CD4 and CD8 B-Gal+ Cells",
  main.cex = 1.5,
  main.pos = c(0.5, 0.95),   # Adjust vertical position of title (0–1)
  main.just = c("center", "top")  # Keep title closer to plot
)

# ---- 7. Save Venn Plot ----
png(file.path(dge_output_dir, "Venn_CD4_CD8_BGal_Pos.png"), width = 1800, height = 1800, res = 300)
grid.draw(venn.plot)
dev.off()

########################## Heatmap ####################################################

sm_genes <- c(
  "S100A9", "CXCL8", "CST3", "TYROBP", "LST1", "FCER1G", "LYZ", "CCL3", "S100A8",
  "CTSS", "AIF1", "S100A12", "SAT1", "G0S2", "S100A11", "PSAP", "NEAT1", "CSTA", "SERPINA1"
)

# ---- 1. Get normalized expression (rlog or counts) ----
rld <- rlog(dds_all, blind = TRUE)
norm_expr <- assay(rld)  # genes x samples

# ---- 2. Filter metadata and add Group label ----
metadata$Group <- paste0(metadata$Cell_Type, "_", metadata$B_Gal)
metadata_rld <- metadata[colnames(norm_expr), ]

# ---- 3. Map Ensembl to gene symbols ----
# Get mapping from previous biomaRt call or redo:
symbol_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(norm_expr),
  mart = ensembl
)

# Match gene symbols to rlog matrix
rownames(norm_expr) <- symbol_map$hgnc_symbol[match(rownames(norm_expr), symbol_map$ensembl_gene_id)]

# ---- 4. Subset for sm_genes ----
norm_subset <- norm_expr[rownames(norm_expr) %in% sm_genes, ]

# ---- 5. Average expression per group ----
metadata_rld <- as_tibble(metadata_rld, rownames = "Sample")

group_expr <- sapply(split(metadata_rld$Sample, metadata_rld$Group), function(samples) {
  rowMeans(norm_subset[, samples, drop = FALSE], na.rm = TRUE)
})

# ---- 6. Create Heatmap ----
pheatmap(
  group_expr,
  cluster_rows = T,
  cluster_cols = T,
  scale = "row",
  color = colorRampPalette(c("forestgreen", "white", "navy"))(100),
  main = "Expression of Selected Genes Across CD4/CD8 B-Gal+/-",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA
)

# Optional: save to file

png(file = file.path(dge_output_dir, "Heatmap_sm_genes_CD4_CD8_BGal.png"), 
    width = 2400, height = 3200, res = 300)

pheatmap(
  group_expr,
  cluster_rows = T,
  cluster_cols = T,
  scale = "row",
  color = colorRampPalette(c("forestgreen", "white", "navy"))(100),
  main = "Expression of Selected Genes Across CD4/CD8 B-Gal+/-",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA
)

dev.off()



#### Heatmap of Shared genes
# shared_genes <- shared_genes[-1]

# ---- 1. Get normalized expression (rlog or counts) ----
rld <- rlog(dds_all, blind = TRUE)
norm_expr <- assay(rld)  # genes x samples

# ---- 2. Filter metadata and add Group label ----
metadata$Group <- paste0(metadata$Cell_Type, "_", metadata$B_Gal)
metadata_rld <- metadata[colnames(norm_expr), ]

# ---- 3. Map Ensembl to gene symbols ----
# Get mapping from previous biomaRt call or redo:
symbol_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(norm_expr),
  mart = ensembl
)

##  Match gene symbols to rlog matrix ----
rownames(norm_expr) <- symbol_map$hgnc_symbol[match(rownames(norm_expr), symbol_map$ensembl_gene_id)]

# ---- 4. Subset for shared_genes ----
norm_subset <- norm_expr[rownames(norm_expr) %in% shared_genes, ]

# ---- 5. Average expression per group ----
metadata_rld <- as_tibble(metadata_rld, rownames = "Sample")

group_expr <- sapply(split(metadata_rld$Sample, metadata_rld$Group), function(samples) {
  rowMeans(norm_subset[, samples, drop = FALSE], na.rm = TRUE)
})

# ---- 6. Create Heatmap ----
pheatmap(
  group_expr,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("forestgreen", "white", "navy"))(100),
  main = "Expression of Shared DEGs Across CD4/CD8 B-Gal+/-",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA
)

# Generate heatmap with custom color
png(file = file.path(dge_output_dir, "Heatmap_SharedGenes_CD4_CD8_BGal.png"),
    width = 2400, height = 3200, res = 300)

pheatmap(
  group_expr,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("forestgreen", "white", "navy"))(100),
  main = "Expression of Shared DEGs Across CD4/CD8 B-Gal+/-",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA
)

dev.off()


################ Heatmaps for top DGE ################################
symbol_map <- symbol_map %>% rename(Gene = ensembl_gene_id, Symbol = hgnc_symbol)
symbol_map <- symbol_map %>% rename(Ensembl_ID = Gene)

# Join to CD4 and CD8 result tables
res_cd4_df <- left_join(res_cd4_df, symbol_map, by = "Ensembl_ID")
res_cd8_df <- left_join(res_cd8_df, symbol_map, by = "Ensembl_ID")

res_cd4_df <- res_cd4_df %>% rename(Gene = Ensembl_ID)
res_cd8_df <- res_cd8_df %>% rename(Gene = Ensembl_ID)

# Annotate CD4 results
res_cd4_df$symbol <- symbol_map$hgnc_symbol[match(res_cd4_df$Gene, symbol_map$ensembl_gene_id)]

# Annotate CD8 results
res_cd8_df$symbol <- symbol_map$hgnc_symbol[match(res_cd8_df$Gene, symbol_map$ensembl_gene_id)]

# ---- 1. Select Top 20 Up and Down Genes for CD4 ----
# Only keep genes with non-blank symbols
res_cd4_df_filtered <- res_cd4_df %>%
  filter(Symbol != "" & !is.na(padj))

# Top 20 upregulated
top_cd4_up <- res_cd4_df_filtered %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

# Top 20 downregulated
top_cd4_down <- res_cd4_df_filtered %>%
  arrange(log2FoldChange) %>%
  head(20)

# Combine gene symbols
cd4_top_genes <- unique(c(top_cd4_up$Symbol, top_cd4_down$Symbol))

# ---- 2. Subset expression matrix ----
norm_cd4 <- norm_expr[rownames(norm_expr) %in% cd4_top_genes, ]
metadata_cd4_rld <- metadata_rld %>% filter(Cell_Type == "CD4")

group_expr_cd4 <- sapply(
  split(metadata_cd4_rld$Sample, metadata_cd4_rld$Group),
  function(samples) rowMeans(norm_cd4[, samples, drop = FALSE], na.rm = TRUE)
)

# ---- 3. Heatmap for CD4 ----

png(file = file.path(dge_output_dir, "Heatmap_TopGenes_CD4_BGal.png"),
    width = 2400, height = 3200, res = 300)

pheatmap(
  group_expr_cd4,
  cluster_rows = T,
  cluster_cols = T,
  color = colorRampPalette(c("forestgreen", "white", "navy"))(100),
  scale='column',
  main = "Top DEGs in CD4 B-Gal+ vs B-Gal-",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA
)

dev.off()

# ---- 4. Repeat for CD8 ----

# ---- 1. Select Top 20 Up and Down Genes for CD4 ----
# Only keep genes with non-blank symbols
res_cd8_df_filtered <- res_cd8_df %>%
  filter(Symbol != "" & !is.na(padj))

# Top 20 upregulated
top_cd8_up <- res_cd4_df_filtered %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

# Top 20 downregulated
top_cd8_down <- res_cd4_df_filtered %>%
  arrange(log2FoldChange) %>%
  head(20)

# Combine gene symbols
cd8_top_genes <- unique(c(top_cd8_up$Symbol, top_cd8_down$Symbol))

# ---- 2. Subset expression matrix ----
norm_cd8 <- norm_expr[rownames(norm_expr) %in% cd8_top_genes, ]
metadata_cd8_rld <- metadata_rld %>% filter(Cell_Type == "CD8")

group_expr_cd8 <- sapply(
  split(metadata_cd8_rld$Sample, metadata_cd8_rld$Group),
  function(samples) rowMeans(norm_cd8[, samples, drop = FALSE], na.rm = TRUE)
)

# ---- 3. Heatmap for CD8 ----

png(file = file.path(dge_output_dir, "Heatmap_TopGenes__CD8_BGal.png"),
    width = 2400, height = 3200, res = 300)

pheatmap(
  group_expr_cd8,
  cluster_rows = T,
  cluster_cols = T,
  color = colorRampPalette(c("forestgreen", "white", "navy"))(100),
  scale='column', ### SCALE BY ROW FOR BLOCK
  main = "Top DEGs in CD8 B-Gal+ vs B-Gal-",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA
)

dev.off()
