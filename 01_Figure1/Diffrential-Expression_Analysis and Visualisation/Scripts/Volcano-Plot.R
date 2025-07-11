# volcano_plot.R
# Generate volcano plots for differential expression analysis
# Use different input files for OSNs or Fatbody as needed

library(readr)
library(ggplot2)
library(scales)
library(dplyr)
library(ggrepel)

# Input file 
# Replace with either "InR_Fatbody_All.csv" or "InR_OSNs_All.csv"
data <- read_csv("InR_Fatbody_All.csv")

# Thresholds 
log2fc_cutoff <- 1
padj_cutoff <- 0.05

#  Categorize regulation 
data$Regulation <- factor(
  ifelse(data$log2FoldChange <= -log2fc_cutoff & data$padj < padj_cutoff, 'Downregulated',
         ifelse(data$log2FoldChange >= log2fc_cutoff & data$padj < padj_cutoff, 'Upregulated', 'Not Significant')),
  levels = c('Upregulated', 'Downregulated', 'Not Significant')
)

# Color scheme 
colors <- c('Upregulated' = '#E41A1C', 'Downregulated' = '#6A9FD8', 'Not Significant' = 'gray70')

#  Filter out NA
data <- data %>% filter(!is.na(log2FoldChange) & !is.na(padj))

#  Select top genes 
top_up <- data %>% filter(Regulation == "Upregulated") %>% arrange(desc(log2FoldChange)) %>% head(10)
top_down <- data %>% filter(Regulation == "Downregulated") %>% arrange(log2FoldChange) %>% head(10)
top_genes <- bind_rows(top_up, top_down)

#  Plot 
plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(-7, 7, by = 1), limits = c(-7, 7)) +
  scale_y_continuous(limits = c(0, 20), oob = squish) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = 'Volcano Plot of Differential Gene Expression',
    subtitle = 'Significance Threshold: |Log2FC| > 1 & padj < 0.05',
    x = expression(Log[2] ~ Fold ~ Change),
    y = expression(-Log[10] ~ Adjusted ~ p-value),
    color = 'Regulation'
  ) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black", size = 0.8) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black", size = 0.8) +
  geom_text_repel(data = top_genes, aes(label = gene_symbol), size = 4, color = "black", box.padding = 0.5)

# Save plot 
ggsave("Volcano_Plot_Fatbody.pdf", plot, width = 8, height = 8)
ggsave("Volcano_Plot_Fatbody.png", plot, width = 8, height = 8, dpi = 300)
