# MSigDB Enrichment Analysis for Drosophila RNA-seq
# Description: Performs MSigDB-based pathway enrichment (KEGG subset) for up/downregulated genes in OSNs and Fatbody datasets

library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Dm.eg.db)
library(msigdbr)
library(stringdist)

# Input and output paths
fat_file <- "C:/Gene_Analysis/InR_Fatbody_All.csv"
osn_file <- "C:/Gene_Analysis/InR_OSNs_All.csv"
output_dir <- "C:/Gene_Analysis/Pathway_Enrichment"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Retrieve Drosophila KEGG gene sets from MSigDB
get_msigdb_kegg_sets <- function() {
  msigdbr(species = "Drosophila melanogaster", category = "C2") %>%
    filter(grepl("KEGG", gs_name, ignore.case = TRUE)) %>%
    select(term = gs_name, gene = entrez_gene) %>%
    mutate(term = as.character(term), gene = as.character(gene))
}

# Map SYMBOL → ENTREZID
map_symbols_to_entrez <- function(symbols) {
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db) %>%
    pull(ENTREZID) %>% unique() %>% as.character()
}

# Map ENTREZID → SYMBOL
entrez_to_symbol <- function(ids) {
  mapping <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Dm.eg.db)
  setNames(mapping$SYMBOL, mapping$ENTREZID)
}

# Annotate gene symbols in enrichment results
add_symbol_column <- function(df) {
  all_entrez <- unique(unlist(strsplit(df$geneID, "/")))
  map <- entrez_to_symbol(all_entrez)
  df$GeneList <- sapply(strsplit(df$geneID, "/"), function(ids) {
    paste(na.omit(map[ids]), collapse = ", ")
  })
  df
}

# Collapse similar pathway terms using Jaccard-like string similarity
filter_similar_terms <- function(df, threshold = 0.8) {
  df <- df[order(df$p.adjust), ]
  keep <- logical(nrow(df))
  used <- rep(FALSE, nrow(df))
  terms <- df$Description
  
  for (i in seq_along(terms)) {
    if (!used[i]) {
      keep[i] <- TRUE
      distances <- stringdist::stringsim(terms[i], terms[(i+1):length(terms)], method = "jw")
      used[(i+1):length(terms)][distances > threshold] <- TRUE
    }
  }
  df[keep, ]
}

# Run enrichment for a dataset (OSNs or Fatbody)
process_dataset <- function(data_file, dataset_name, msigdb_sets) {
  data <- read_csv(data_file, show_col_types = FALSE)
  
  up_genes <- data %>% filter(padj < 0.05 & log2FoldChange > 0) %>% pull(gene_symbol) %>% unique()
  down_genes <- data %>% filter(padj < 0.05 & log2FoldChange < 0) %>% pull(gene_symbol) %>% unique()
  
  all_results <- list()
  
  for (condition in c("up", "down")) {
    gene_symbols <- if (condition == "up") up_genes else down_genes
    entrez_ids <- map_symbols_to_entrez(gene_symbols)
    
    if (length(entrez_ids) > 0) {
      enrichment <- enricher(gene = entrez_ids, TERM2GENE = msigdb_sets, pvalueCutoff = 0.05)
      
      if (!is.null(enrichment) && nrow(enrichment@result) > 0) {
        df <- add_symbol_column(as.data.frame(enrichment))
        df <- df %>% mutate(Dataset = dataset_name, Regulation = condition)
        write_csv(df, file.path(output_dir, paste0("MSigDB_", dataset_name, "_", condition, ".csv")))
        all_results[[condition]] <- df
      }
    }
  }
  
  if (length(all_results) > 0) {
    filtered <- bind_rows(all_results) %>%
      group_by(Dataset, Regulation) %>%
      group_modify(~ filter_similar_terms(.x))
    write_csv(filtered, file.path(output_dir, paste0("MSigDB_", dataset_name, "_Filtered.csv")))
  }
}

# Run enrichment for both OSNs and Fatbody
msigdb_sets <- get_msigdb_kegg_sets()
process_dataset(fat_file, "Fatbody", msigdb_sets)
process_dataset(osn_file, "OSNs", msigdb_sets)
