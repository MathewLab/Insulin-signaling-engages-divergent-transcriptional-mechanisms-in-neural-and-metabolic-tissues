# GO Enrichment Analysis using topGO with elim Fisher test
# For OSNs and Fatbody insulin knockdown RNA-seq data
# Outputs: CSVs with adjusted GO terms (padj < 0.05)

# Load required libraries
library(readr)
library(dplyr)
library(topGO)
library(org.Dm.eg.db)


# Input data (relative path)
fat_file <- "data/InR_Fatbody_All.csv"
osn_file <- "data/InR_OSNs_All.csv"
output_dir <- "results/"
if (!dir.exists(output_dir)) dir.create(output_dir)


# Load differential expression data
fatbody <- read_csv(fat_file, show_col_types = FALSE)
osns <- read_csv(osn_file, show_col_types = FALSE)

# Define significant gene sets
fat_genes_sig <- fatbody %>%
  filter(!is.na(padj), padj < 0.05) %>%
  pull(gene_symbol) %>% unique()

osns_genes_sig <- osns %>%
  filter(!is.na(padj), padj < 0.05) %>%
  pull(gene_symbol) %>% unique()

# Background gene sets
all_fat_genes <- unique(fatbody$gene_symbol)
all_osns_genes <- unique(osns$gene_symbol)

# Enrichment function
run_topgo_analysis <- function(gene_list, all_genes, dataset_name) {
  gene_universe <- factor(as.integer(all_genes %in% gene_list))
  names(gene_universe) <- all_genes

  go_data <- new("topGOdata",
                 ontology = "BP",
                 allGenes = gene_universe,
                 geneSelectionFun = function(x) x == 1,
                 annot = annFUN.org,
                 mapping = "org.Dm.eg.db",
                 ID = "SYMBOL")

  result_elim <- runTest(go_data, algorithm = "elim", statistic = "fisher")
  all_res <- score(result_elim)
  go_terms <- names(all_res)

  # Adjust p-values
  pval_df <- data.frame(
    GO.ID = go_terms,
    raw.pval = all_res,
    adj.pval = p.adjust(all_res, method = "BH")
  )

  # Add term descriptions
  pval_df$Term <- Term(GOTERM[go_terms])

  # Filter significant terms
  filtered <- pval_df %>%
    filter(adj.pval < 0.05) %>%
    arrange(adj.pval) %>%
    select(GO.ID, Term, raw.pval, adj.pval)

  # Write output
  output_file <- file.path(output_dir, paste0("GO_Results_", dataset_name, "_elim_fisher_adj.csv"))
  write_csv(filtered, output_file)

  return(filtered)
}

# Run enrichment for both datasets
go_fat <- run_topgo_analysis(fat_genes_sig, all_fat_genes, "Fatbody")
go_osns <- run_topgo_analysis(osns_genes_sig, all_osns_genes, "OSNs")
