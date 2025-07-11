library(readr) 
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Dm.eg.db)
library(stringdist)

# Define input and output
fat_file <- "C:/Gene_Analysis/InR_Fatbody_All.csv"
osn_file <- "C:/Gene_Analysis/InR_OSNs_All.csv"
output_dir <- "C:/Gene_Analysis/Pathway_Enrichment"
final_summary_file <- file.path(output_dir, "Final_Filtered_Pathways.csv")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Map SYMBOL to ENTREZID and KEGG
get_entrez_ids <- function(symbols) {
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
}

# Map ENTREZID to SYMBOL
entrez_to_symbol <- function(ids) {
  map <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Dm.eg.db)
  setNames(map$SYMBOL, map$ENTREZID)
}

# Enrichment wrapper
run_enrichment <- function(entrez_ids, type = "KEGG") {
  if (type == "KEGG") {
    enrichKEGG(gene = entrez_ids, organism = 'dme', keyType = "ncbi-geneid", pvalueCutoff = 0.05)
  } else if (type == "Reactome") {
    enrichPathway(gene = entrez_ids, organism = "fly", pvalueCutoff = 0.05)
  }
}

# Add gene list column (SYMBOL version)
add_gene_list_column <- function(enrich_res) {
  enrich_df <- as.data.frame(enrich_res)
  entrez_ids <- unique(unlist(strsplit(enrich_df$geneID, "/")))
  symbol_map <- entrez_to_symbol(entrez_ids)
  enrich_df$GeneList <- sapply(strsplit(enrich_df$geneID, "/"), function(ids) {
    paste(na.omit(symbol_map[ids]), collapse = ", ")
  })
  enrich_df
}

# Filter similar pathway terms using Jaccard-like distance
filter_similar_pathways <- function(df, similarity_threshold = 0.8) {
  keep <- logical(nrow(df))
  df <- df[order(df$p.adjust), ]
  terms <- df$Description
  used <- rep(FALSE, length(terms))
  
  for (i in seq_along(terms)) {
    if (!used[i]) {
      keep[i] <- TRUE
      distances <- stringdist::stringsim(terms[i], terms[(i+1):length(terms)], method = "jw")
      used[(i+1):length(terms)][distances > similarity_threshold] <- TRUE
    }
  }
  
  df[keep, ]
}

# Dataset processor
process_dataset <- function(file, dataset) {
  data <- read_csv(file, show_col_types = FALSE)
  
  # Filter genes
  up_genes <- data %>% filter(!is.na(padj) & padj < 0.05 & log2FoldChange > 0) %>% pull(gene_symbol) %>% unique()
  down_genes <- data %>% filter(!is.na(padj) & padj < 0.05 & log2FoldChange < 0) %>% pull(gene_symbol) %>% unique()
  
  # Convert to ENTREZ
  up_entrez <- get_entrez_ids(up_genes)
  down_entrez <- get_entrez_ids(down_genes)
  
  final_combined <- list()
  
  for (pathway_type in c("KEGG", "Reactome")) {
    for (group in c("up", "down")) {
      genes <- if (group == "up") up_entrez else down_entrez
      
      if (length(genes) > 0) {
        enrich_res <- tryCatch({
          run_enrichment(genes, type = pathway_type)
        }, error = function(e) NULL)
        
        if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
          enriched_df <- add_gene_list_column(enrich_res)
          enriched_df <- enriched_df %>% mutate(Dataset = dataset, Regulation = group, Source = pathway_type)
          
          fname <- paste0(pathway_type, "_", dataset, "_", group, ".csv")
          write.csv(enriched_df, file.path(output_dir, fname), row.names = FALSE)
          
          final_combined[[paste(dataset, group, pathway_type)]] <- enriched_df
        }
      }
    }
  }
  
  # Combine all results and filter similar terms
  combined_df <- bind_rows(final_combined)
  if (nrow(combined_df) > 0) {
    filtered <- combined_df %>% group_by(Dataset, Regulation, Source) %>% group_modify(~ filter_similar_pathways(.x))
    write.csv(filtered, final_summary_file, row.names = FALSE)
  }
}

# Run enrichment
process_dataset(fat_file, "Fatbody")
process_dataset(osn_file, "OSNs")
