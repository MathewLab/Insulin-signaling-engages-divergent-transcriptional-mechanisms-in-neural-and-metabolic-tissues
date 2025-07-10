# venn_diagram.R
# Create Venn diagrams comparing DEGs between OSNs and Fatbody
library(readr)
library(dplyr)
library(VennDiagram)
library(grid)

# Input files 
fatbody_data <- read_csv("InR_Fatbody_All.csv")
osns_data <- read_csv("InR_OSNs_All.csv")

#  Venn plot function
generate_venn <- function(set1, set2, filename, title) {
  venn.plot <- draw.pairwise.venn(
    area1 = length(set1),
    area2 = length(set2),
    cross.area = length(intersect(set1, set2)),
    category = c("Fatbody", "OSNs"),
    fill = c("#66B2FF", "#FF6666"),
    alpha = 0.6,
    lwd = 2,
    cex = 2,
    cat.cex = 2,
    cat.fontface = "bold",
    cat.pos = c(-30, 30),
    cat.dist = 0.05
  )

  png(paste0(filename, ".png"), width = 1200, height = 1200, res = 300)
  grid.draw(venn.plot)
  dev.off()

  pdf(paste0(filename, ".pdf"), width = 7, height = 7)
  grid.draw(venn.plot)
  dev.off()
}

#  Extract significant gene sets 
fatbody_padj_0.05 <- fatbody_data %>% filter(padj < 0.05) %>% pull(gene_symbol)
osns_padj_0.05 <- osns_data %>% filter(padj < 0.05) %>% pull(gene_symbol)

fatbody_padj_0.01 <- fatbody_data %>% filter(padj < 0.01) %>% pull(gene_symbol)
osns_padj_0.01 <- osns_data %>% filter(padj < 0.01) %>% pull(gene_symbol)

fatbody_fc_padj <- fatbody_data %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_symbol)
osns_fc_padj <- osns_data %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_symbol)

# Generate plots 
generate_venn(fatbody_padj_0.05, osns_padj_0.05, "Venn_padj_0.05", "padj < 0.05")
generate_venn(fatbody_padj_0.01, osns_padj_0.01, "Venn_padj_0.01", "padj < 0.01")
generate_venn(fatbody_fc_padj, osns_fc_padj, "Venn_padj_FC", "padj < 0.05 & |log2FC| > 1")
