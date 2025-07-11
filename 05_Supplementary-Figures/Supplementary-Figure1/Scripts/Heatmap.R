# Heatmap of Top 30 Shared DEGs per Regulatory Pattern (Vertical View)
# Description: Generates quadrant-specific vertical heatmaps (log2FC) for shared DEGs in OSNs and Fat Bod

library(readr)
library(dplyr)
library(tibble)
library(pheatmap)

# Custom diverging color palette
my_colors <- colorRampPalette(c("#0066CC", "white", "#CC3300"))(100)

# Load input data
osn <- read_csv("InR_OSNs_All.csv")
fat <- read_csv("InR_Fatbody_All.csv")

#  Filter significant DEGs (padj < 0.05 and |log2FC| > 1)
osn_sig <- osn %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  select(gene_symbol, log2FoldChange, padj) %>%
  rename(log2FC_OSNs = log2FoldChange, padj_OSNs = padj)

fat_sig <- fat %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  select(gene_symbol, log2FoldChange, padj) %>%
  rename(log2FC_Fatbody = log2FoldChange, padj_Fatbody = padj)

#  Merge and annotate regulation patterns
common_genes <- inner_join(osn_sig, fat_sig, by = "gene_symbol") %>%
  mutate(
    regulation_pattern = case_when(
      log2FC_OSNs > 0 & log2FC_Fatbody > 0 ~ "Up.Up",
      log2FC_OSNs < 0 & log2FC_Fatbody < 0 ~ "Down.Down",
      log2FC_OSNs > 0 & log2FC_Fatbody < 0 ~ "Up.Down",
      log2FC_OSNs < 0 & log2FC_Fatbody > 0 ~ "Down.Up"
    )
  ) %>%
  filter(!is.na(regulation_pattern))

# Regulation patterns
patterns <- c("Up.Up", "Down.Down", "Up.Down", "Down.Up")

# Select top 30 genes per pattern based on combined padj score
get_top30_by_padj <- function(data, pattern) {
  data %>%
    filter(regulation_pattern == pattern) %>%
    mutate(score = -log10(padj_OSNs + padj_Fatbody)) %>%
    slice_max(order_by = score, n = 30)
}

# Create log2FC matrix for pheatmap
prepare_log2fc_vertical <- function(df) {
  df %>%
    mutate(gene_label = gene_symbol) %>%
    select(gene_label, log2FC_OSNs, log2FC_Fatbody) %>%
    distinct() %>%
    column_to_rownames("gene_label") %>%
    as.matrix()
}

# Color scale settings (shared across all heatmaps)
fc_breaks <- seq(-4, 4, length.out = 101)

# Save PDF and PNG for each vertical heatmap
save_fc_vertical <- function(mat, pattern) {
  safe_pattern <- gsub("\\.", "_", pattern)
  filename_base <- paste0("Heatmap_log2FC_", safe_pattern, "_Vertical")
  
  # PDF output
  pheatmap(mat,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           color = my_colors,
           breaks = fc_breaks,
           border_color = NA,
           angle_col = 45,
           show_colnames = TRUE,
           fontsize = 8,
           cellwidth = 12,
           cellheight = 12,
           main = paste("log2FC Vertical Heatmap:", pattern),
           filename = paste0(filename_base, ".pdf"),
           width = 4, height = 7)

  # PNG output
  png(paste0(filename_base, ".png"), width = 400, height = 700, res = 100)
  pheatmap(mat,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           color = my_colors,
           breaks = fc_breaks,
           border_color = NA,
           angle_col = 45,
           show_colnames = TRUE,
           fontsize = 8,
           cellwidth = 12,
           cellheight = 12,
           main = paste("log2FC Vertical Heatmap:", pattern),
           filename = NA)
  dev.off()
}

# Generate and save heatmaps for each regulation pattern
for (p in patterns) {
  top_genes <- get_top30_by_padj(common_genes, p)
  mat_fc <- prepare_log2fc_vertical(top_genes)
  save_fc_vertical(mat_fc, p)
}
