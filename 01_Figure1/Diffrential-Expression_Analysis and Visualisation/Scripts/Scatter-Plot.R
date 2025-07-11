# Scatter-Plot.R
# Scatter-quadrant-Plot of log2 fold changes for shared DEGs for directionality (padj < 0.05)

library(readr)
library(dplyr)
library(ggplot2)
library(scales)

#  Input 
osn <- read_csv("InR_OSNs_All.csv")
fat <- read_csv("InR_Fatbody_All.csv")

# Filter padj < 0.05 
osn_sig <- osn %>% filter(padj < 0.05) %>%
  select(gene_symbol, log2FoldChange, padj) %>%
  rename(log2FC_OSNs = log2FoldChange, padj_OSNs = padj)

fat_sig <- fat %>% filter(padj < 0.05) %>%
  select(gene_symbol, log2FoldChange, padj) %>%
  rename(log2FC_Fatbody = log2FoldChange, padj_Fatbody = padj)

# Join datasets 
common_genes <- inner_join(osn_sig, fat_sig, by = "gene_symbol") %>%
  mutate(
    reg_OSNs = ifelse(log2FC_OSNs > 0, "Up", "Down"),
    reg_Fatbody = ifelse(log2FC_Fatbody > 0, "Up", "Down"),
    regulation_pattern = paste0(reg_OSNs, ".", reg_Fatbody)
  )

# Save filtered genes 
write_csv(common_genes, "Common_Genes_padj0.05_OSNs_Fatbody.csv")

# Pearson correlation
cor_result <- cor.test(common_genes$log2FC_OSNs, common_genes$log2FC_Fatbody)
cor_label <- paste0("r = ", round(cor_result$estimate, 2),
                    ", p = ", format.pval(cor_result$p.value, digits = 2, scientific = TRUE))

# Plot 
pairwise_plot <- ggplot(common_genes, aes(x = log2FC_OSNs, y = log2FC_Fatbody)) +
  geom_point(aes(color = regulation_pattern), alpha = 0.7, size = 2.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  annotate("text", x = 1.5, y = 2.7, label = cor_label, size = 5, hjust = 0, color = "black") +
  scale_x_continuous(limits = c(-4, 4), oob = squish) +
  scale_y_continuous(limits = c(-4, 4), oob = squish) +
  annotate("rect", xmin = -3, xmax = 3, ymin = -3, ymax = 3, color = "gray80", fill = NA, linetype = "dashed") +
  scale_color_manual(values = c("Up.Up" = "red2", "Down.Down" = "#1E90FF", "Up.Down" = "orange", "Down.Up" = "purple")) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(), plot.margin = margin(30, 30, 30, 30)) +
  labs(title = "Regulation Directionality of Common Genes (padj < 0.05)",
       x = "log2 Fold Change (OSNs)", y = "log2 Fold Change (Fatbody)")

ggsave("Regulation_Directionality_of_CommonGenes.pdf", pairwise_plot, width = 7.5, height = 7.5, dpi = 300)
