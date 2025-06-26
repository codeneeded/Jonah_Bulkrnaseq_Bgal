library(dplyr)
library(enrichR)
library(openxlsx)
library(ggplot2)

# Set Enrichr databases
databases <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019"
)

tf_databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- setdiff(databases, tf_databases)

# Output directory
base_output <- "/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/Enrichr/"
dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# Read CSVs explicitly
cd4_dge <- read.csv("/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/DGE/CD4_BGal_Pos_vs_Neg_DGE.csv")
cd8_dge <- read.csv("/home/akshay-iyer/Documents/Jonah_Bulkrnaseq_Bgal/DGE/CD8_BGal_Pos_vs_Neg_DGE.csv")

run_enrichment <- function(dge_df, label, direction) {
  if (nrow(dge_df) == 0) {
    message("No genes found for ", label, " (", direction, ")")
    return(NULL)
  }
  
  gene_list <- unique(dge_df$hgnc_symbol[dge_df$hgnc_symbol != ""])
  
  if (length(gene_list) == 0) {
    message("No valid gene symbols for ", label, " (", direction, ")")
    return(NULL)
  }
  
  enrichment <- enrichr(gene_list, databases)
  
  # Save Excel output
  excel_out <- file.path(base_output, paste0(label, "_", direction, "_Enrichment.xlsx"))
  wb <- createWorkbook()
  for (db in names(enrichment)) {
    addWorksheet(wb, substr(db, 1, 31))
    writeData(wb, substr(db, 1, 31), enrichment[[db]])
  }
  saveWorkbook(wb, excel_out, overwrite = TRUE)
  
  # Collect significant terms
  top_tf_list <- list()
  top_pathway_list <- list()
  
  for (db_name in names(enrichment)) {
    db_results <- enrichment[[db_name]]
    
    # Rename if necessary
    if ("Combined Score" %in% colnames(db_results)) {
      db_results <- db_results %>% rename(Combined.Score = `Combined Score`)
    }
    
    if (!"Combined.Score" %in% colnames(db_results)) {
      message("Skipping ", db_name, " for ", label, " (", direction, ") — no Combined.Score.")
      next
    }
    
    sig_results <- db_results %>% filter(Adjusted.P.value < 0.05)
    
    if (nrow(sig_results) > 0) {
      top_terms <- sig_results %>%
        arrange(desc(Combined.Score)) %>%
        slice_head(n = 10) %>%
        mutate(Database = db_name)
      
      if (db_name %in% tf_databases) {
        top_tf_list[[db_name]] <- top_terms
      } else if (db_name %in% pathway_databases) {
        top_pathway_list[[db_name]] <- top_terms
      }
    }
  }
  
  # Process TFs
  tf_df <- bind_rows(top_tf_list)
  if ("Combined.Score" %in% colnames(tf_df) && nrow(tf_df) > 0) {
    tf_df <- tf_df %>%
      arrange(desc(Combined.Score)) %>%
      slice_head(n = 20)
    
    p_tf <- ggplot(tf_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        title = paste("Top Transcription Factors -", label, direction),
        x = "TF Term", y = "Combined Score"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    
    ggsave(
      filename = file.path(base_output, paste0(label, "_", direction, "_Transcription_Factors.png")),
      plot = p_tf, width = 12, height = 10, dpi = 300, bg = "white"
    )
  } else {
    message("No valid TF results for ", label, " (", direction, ")")
  }
  
  # Process pathways
  pathway_df <- bind_rows(top_pathway_list)
  if ("Combined.Score" %in% colnames(pathway_df) && nrow(pathway_df) > 0) {
    pathway_df <- pathway_df %>%
      arrange(desc(Combined.Score)) %>%
      slice_head(n = 20)
    
    p_path <- ggplot(pathway_df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = Database)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        title = paste("Top Pathways -", label, direction),
        x = "Pathway Term", y = "Combined Score"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    
    ggsave(
      filename = file.path(base_output, paste0(label, "_", direction, "_Pathways.png")),
      plot = p_path, width = 12, height = 10, dpi = 300, bg = "white"
    )
  } else {
    message("No valid pathway results for ", label, " (", direction, ")")
  }
}


# Define filtering and call function for each case
process_comparison <- function(dge_df, label) {
  up_df   <- dge_df %>% filter(padj < 0.05 & log2FoldChange > 1)
  down_df <- dge_df %>% filter(padj < 0.05 & log2FoldChange < -1)
  all_df  <- dge_df %>% filter(padj < 0.05)
  
  run_enrichment(up_df,   label, "Up")
  run_enrichment(down_df, label, "Down")
  run_enrichment(all_df,  label, "All")
}

# Run for both CD4 and CD8
process_comparison(cd4_dge, "CD4_BGal_Pos_vs_Neg")
process_comparison(cd8_dge, "CD8_BGal_Pos_vs_Neg")

