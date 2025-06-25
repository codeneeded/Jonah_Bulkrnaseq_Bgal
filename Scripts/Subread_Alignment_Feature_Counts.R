#### Libraries
library(Rsubread)
library(readxl)
library(writexl)
library(openxlsx)

############################################################## Building Index File ############################################################

# Path to FASTA file (unzipped)
fasta <- "~/Documents/Jonah_Bulkrnaseq_Bgal/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
index_basename <- "~/Documents/Jonah_Bulkrnaseq_Bgal/Index_hg38/hg38_index"

# Build index (not split — full index)
# buildindex(basename = index_basename, reference = fasta, gappedIndex = FALSE, indexSplit = FALSE)

################################################################### Alignment + Feature Counts ############################################


samples <- sprintf("kupritz-26037-%03d", 1:12) # Define sample ID


### Read in Path Names

fastq_dir <- "~/Documents/Jonah_Bulkrnaseq_Bgal/FastQ_Trimmed"
index_basename <- "~/Documents/Jonah_Bulkrnaseq_Bgal/Index_hg38/hg38_index"
output_dir <- "~/Documents/Jonah_Bulkrnaseq_Bgal/Feature Counts"
gtf_file <- "~/Documents/Jonah_Bulkrnaseq_Bgal/Homo_sapiens.GRCh38.114.gtf"


#### Alignment and Feature Count- For LOOP

if (!dir.exists(output_dir)) dir.create(output_dir) ## Create output folder if it doesnt exist

# Initialize list to collect counts
all_counts_list <- list() # Empty List to store output

for (sample in samples) {
  cat("\n🔄 Processing sample:", sample, "\n")
  
  r1 <- file.path(fastq_dir, paste0(sample, "_clipr_r1.fastq.gz"))
  r2 <- file.path(fastq_dir, paste0(sample, "_clipr_r2.fastq.gz"))
  bam_file <- file.path(output_dir, paste0(sample, ".bam"))
  
  align(
    index = index_basename,
    readfile1 = r1,
    readfile2 = r2,
    input_format = "gzFASTQ",
    output_file = bam_file,
    nthreads = 30
  )
  
  fc_results <- featureCounts(
    files = bam_file,
    isPairedEnd = TRUE,
    useMetaFeatures = TRUE,
    annot.ext = gtf_file,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    nthreads = 30
  )
  
  # Store counts with proper column name
  sample_counts <- as.data.frame(fc_results$counts)
  colnames(sample_counts) <- sample
  all_counts_list[[sample]] <- sample_counts
  
  # Save individual sample outputs as before
  write.csv(sample_counts, file = file.path(output_dir, paste0(sample, "_raw_counts.csv")), quote = FALSE)
  write.xlsx(
    list(
      counts = sample_counts,
      annotation = as.data.frame(fc_results$annotation),
      alignment_summary = as.data.frame(fc_results$stat)
    ),
    file = file.path(output_dir, paste0(sample, "_feature_counts.xlsx"))
  )
  
  cat("✔ Done:", sample, "\n")
}

# After the loop — combine all sample counts
combined_counts <- do.call(cbind, all_counts_list)

# Save the combined matrix
write.csv(combined_counts, file = file.path(output_dir, "raw_combined_featureCounts.csv"), quote = FALSE)

