# Bulk RNA-seq of SA-β-Gal⁺ vs SA-β-Gal⁻ Human CD4⁺ and CD8⁺ T Cells

Bulk RNA-seq analysis of **senescence-associated β-galactosidase positive (SA-β-Gal⁺)** versus **negative** human **CD4⁺** and **CD8⁺** T cells. The goal is to characterize the transcriptional program of T-cell senescence, identify shared vs lineage-specific senescence signatures across CD4 and CD8 compartments, and project the resulting signatures onto an external single-cell RNA-seq dataset via module scoring.

SA-β-Gal is the canonical histochemical / flow-cytometric marker of cellular senescence; sorting on SA-β-Gal activity *(e.g. C₁₂FDG or SPiDER-βGal — fill in)* enables transcriptomic comparison of senescent and non-senescent cells from the same donor and lineage.

---

## Experimental design

| Factor | Levels |
|---|---|
| **Cell type** | CD4⁺ T cells, CD8⁺ T cells |
| **SA-β-Gal status** | β-Gal⁺ (senescent), β-Gal⁻ (non-senescent) |
| **Donors** | *(n donors, demographics — fill in)* |
| **Total libraries** | 12 (`kupritz-26037-001` … `-012`) |
| **Source** | *(e.g. PBMCs from healthy donors; sort scheme — fill in)* |
| **Reference** | *Homo sapiens*, GRCh38 *(confirm)* — KEGG hsa codes used downstream |
| **Library prep / sequencing** | *(kit, platform, read length — fill in)* |

Sample-level metadata lives in [`Bulk_RNA-Seq_Metadata.csv`](Bulk_RNA-Seq_Metadata.csv).

---

## Repository structure

```
Jonah_Bulkrnaseq_Bgal/
├── Bulk_RNA-Seq_Metadata.csv             # Sample → cell type / β-Gal status / donor mapping
│
├── Scripts/                              # End-to-end R analysis pipeline
│   ├── Subread_Alignment_Feature_Counts.R   # Rsubread alignment + featureCounts
│   ├── Differential_Expression_Deseq2.R     # DESeq2 DGE (β-Gal⁺ vs β-Gal⁻ within CD4 / CD8)
│   ├── Pathway_Analysis.R                   # enrichR over Up / Down / All DEG sets
│   ├── Kegg_Pathview_Overlay.R              # log2FC overlay on KEGG pathway diagrams
│   └── CD8_Module_Scoring_scRNAseq.R        # Project bulk signatures onto external scRNA-seq
│
├── QC Reports/                           # Per-sample read-level QC (HTML + JSON, kupritz-26037-001…012)
│
├── Feature Counts/                       # Per-sample alignment outputs
│   ├── kupritz-26037-XXX.bam.summary         # featureCounts summary
│   ├── kupritz-26037-XXX.bam.indel.vcf       # Subread indel calls
│   ├── kupritz-26037-XXX_raw_counts.csv      # Raw gene counts
│   ├── kupritz-26037-XXX_feature_counts.xlsx
│   └── raw_combined_featureCounts.csv        # Merged 12-sample count matrix
│
├── DGE/                                  # DESeq2 results & visualizations
│   ├── CD4_BGal_Pos_vs_Neg_DGE.csv
│   ├── CD8_BGal_Pos_vs_Neg_DGE.csv
│   ├── Shared_DEGs_CD4_CD8_BGal_Pos.csv      # DEGs shared between CD4 and CD8 senescence
│   ├── PCA_CellType_BGal.png
│   ├── CD4_BGal_Volcano.png, CD8_BGal_Volcano.png
│   ├── Venn_CD4_CD8_BGal_Pos.png             # Overlap of CD4 / CD8 senescence DEGs
│   ├── Heatmap_TopGenes_CD4_BGal.png, Heatmap_TopGenes_CD8_BGal.png
│   ├── Heatmap_SharedGenes_CD4_CD8_BGal.png  # Shared signature heatmap
│   └── Heatmap_sm_genes_CD4_CD8_BGal.png     # Targeted small-gene-set heatmap
│
├── Enrichr/                              # Functional enrichment per cell type and direction
│   ├── {CD4,CD8}_BGal_Pos_vs_Neg_{Up,Down,All}_Enrichment.xlsx
│   ├── {CD4,CD8}_BGal_Pos_vs_Neg_{Up,Down,All}_Pathways.png
│   └── {CD4,CD8}_BGal_Pos_vs_Neg_{Up,Down,All}_Transcription_Factors.png
│
├── Pathway_Overlay/                      # KEGG pathway diagrams with per-gene log2FC overlay
│   ├── CD4/Apoptosis/hsa04210.CD4_BGal_Pos_vs_Neg.png
│   └── CD8/Apoptosis/hsa04210.CD8_BGal_Pos_vs_Neg.png
│
├── scRNAseq/Module_Score/                # Bulk-derived senescence signatures applied to scRNA-seq
│   ├── ModuleScores_byCluster_sampleMean_A.png
│   └── ModuleScores_byCluster_sampleMean_E.png
│
├── Jonah_Bulkrnaseq_Bgal.Rproj
├── README.md
└── .gitignore
```

---

## Pipeline overview

### 1. Alignment & quantification — `Subread_Alignment_Feature_Counts.R`
FASTQs aligned to GRCh38 with **Rsubread** `align()`; gene-level counts produced with **`featureCounts`** against an Ensembl/GENCODE GTF *(specify version)*. Per-sample summaries, indel VCFs, and raw count CSVs land in `Feature Counts/`. The 12 per-sample matrices are merged into `raw_combined_featureCounts.csv`.

### 2. Read QC — `QC Reports/`
Per-sample HTML + JSON reports (one per library) covering adapter trimming, read-length distributions, and basic alignment QC.

### 3. Differential expression — `Differential_Expression_Deseq2.R`
DESeq2 contrasts run **within each cell type**:

```
CD4: β-Gal⁺ vs β-Gal⁻
CD8: β-Gal⁺ vs β-Gal⁻
```

Outputs: full DEG tables, volcano plots, top-gene heatmaps, a CD4/CD8 PCA colored by β-Gal status, plus a **shared-DEG analysis** (Venn + intersection table + heatmap) to isolate a **core senescence signature** common to both lineages, separated from CD4- or CD8-restricted programs.

### 4. Functional enrichment — `Pathway_Analysis.R`
**enrichR** queried separately for **Up**, **Down**, and **All** DEG lists in each cell type, against pathway and transcription-factor libraries. Excel exports + pathway / TF bar plots in `Enrichr/`. Splitting by direction lets you distinguish what senescent T cells *gain* (e.g. SASP, NF-κB, p53 targets) from what they *lose* (e.g. proliferation programs, ribosome biogenesis).

### 5. KEGG pathway overlay — `Kegg_Pathview_Overlay.R`
**`pathview`** colors KEGG pathway diagrams (organism `hsa`) by per-gene log2FC. Currently focused on **Apoptosis (hsa04210)** for CD4 and CD8 — useful given senescent cells' known apoptosis-resistance phenotype. Add further pathway IDs to the script to extend.

### 6. Single-cell projection — `CD8_Module_Scoring_scRNAseq.R`
The bulk-derived senescence DEG sets are taken into a **scRNA-seq dataset** *(specify which — e.g. an external CD8 atlas)* and applied as **`AddModuleScore`** signatures (Seurat). Per-cluster, sample-mean module scores are written to `scRNAseq/Module_Score/`, providing a single-cell readout of which clusters / cell states most resemble bulk-defined senescent T cells.

---

## Dependencies

**Alignment & counts**
- [`Rsubread`](https://bioconductor.org/packages/Rsubread/) — alignment and `featureCounts`

**DGE & annotation**
- [`DESeq2`](https://bioconductor.org/packages/DESeq2/) — differential expression
- [`org.Hs.eg.db`](https://bioconductor.org/packages/org.Hs.eg.db/), [`AnnotationDbi`](https://bioconductor.org/packages/AnnotationDbi/) — human gene annotation
- `dplyr`, `tidyr`, `readxl`, `openxlsx` — wrangling / IO

**Enrichment & pathway viz**
- [`enrichR`](https://CRAN.R-project.org/package=enrichR) — pathway / TF enrichment
- [`pathview`](https://bioconductor.org/packages/pathview/) — KEGG diagram overlays
- [`clusterProfiler`](https://yulab-smu.top/biomedical-knowledge-mining-book/) *(if used downstream)*

**Visualization**
- [`ComplexHeatmap`](https://bioconductor.org/packages/ComplexHeatmap/), `circlize` — heatmaps
- [`EnhancedVolcano`](https://bioconductor.org/packages/EnhancedVolcano/) — volcano plots
- [`VennDiagram`](https://CRAN.R-project.org/package=VennDiagram) / `ggvenn` — DEG overlap
- `ggplot2`, `patchwork`

**Single-cell projection**
- [`Seurat`](https://satijalab.org/seurat/) — `AddModuleScore` for scRNA-seq projection

---

## Reproducing the analysis

1. Clone and open `Jonah_Bulkrnaseq_Bgal.Rproj` in RStudio.
2. Place raw FASTQs at the path expected by `Scripts/Subread_Alignment_Feature_Counts.R` *(top-of-script `setwd` / `load.path` may need updating to your machine)*.
3. Run scripts in order:
   1. `Subread_Alignment_Feature_Counts.R`
   2. `Differential_Expression_Deseq2.R`
   3. `Pathway_Analysis.R`
   4. `Kegg_Pathview_Overlay.R`
   5. `CD8_Module_Scoring_scRNAseq.R` — requires the external scRNA-seq object *(specify path / accession)*
4. Outputs regenerate under `DGE/`, `Enrichr/`, `Pathway_Overlay/`, `scRNAseq/Module_Score/`.

---

## Key results at a glance

| Question | File(s) |
|---|---|
| What changes in CD4 senescence? | `DGE/CD4_BGal_Pos_vs_Neg_DGE.csv`, `CD4_BGal_Volcano.png` |
| What changes in CD8 senescence? | `DGE/CD8_BGal_Pos_vs_Neg_DGE.csv`, `CD8_BGal_Volcano.png` |
| Core (shared) senescence signature | `DGE/Shared_DEGs_CD4_CD8_BGal_Pos.csv`, `Venn_CD4_CD8_BGal_Pos.png`, `Heatmap_SharedGenes_CD4_CD8_BGal.png` |
| Pathways / TFs driving the program | `Enrichr/*` |
| Apoptosis machinery in senescent T cells | `Pathway_Overlay/{CD4,CD8}/Apoptosis/hsa04210.*` |
| Single-cell correlate of bulk signatures | `scRNAseq/Module_Score/*` |

---
