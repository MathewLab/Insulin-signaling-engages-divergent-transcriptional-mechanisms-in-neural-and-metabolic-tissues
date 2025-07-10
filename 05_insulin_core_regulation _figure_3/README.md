# 04_insulin_core_regulation

This module analyzes the differential impact of reduced insulin signaling on core insulin pathway components in OSNs and Fatbody tissues.

## Overview

This figure combines four analytical components:
1. Venn Diagram of significantly regulated insulin core genes (common vs unique).
2. Pairwise Scatter Plot showing directional regulation (log2FC OSNs vs Fatbody).
3. GO Enrichment analysis and Visualisation using `GO_enrichment_circlize.R`.
4. Reactome-style Pathway Schematic highlighting affected genes and regulation direction.

### 1. Venn Diagram – Insulin Core Components (Figure 3A)
This panel uses the same R script as Figure 1 (03_figure1/scripts/Venn_diagram.R),
but with a filtered list of insulin core genes.

Script reused: ../03_figure1/scripts/Venn_diagram.R Input: data/insulin_pathway_core_gene_list.csv
Output: plots/Venn_InsulinCore.pdf

### 2. Pairwise Scatter Plot – Insulin Core Components (Figure 3B)

This plot shows the direction and strength of log2 fold changes for insulin core genes in OSNs vs Fatbody.
Script reused: `../03_figure1/scripts/pairwise-correlation.R` Input: `data/insulin_pathway_core_gene_list.csv
Output: `plots/Scatter_InsulinCore.pdf`

### 3. GO Enrichment – Circos Visualization (Figure 3C)

GO enrichment analysis was performed using **topGO (v2.56.0)** with the **elim algorithm** and **Fisher’s exact test**.  
Each tissue used its own expressed gene background (≥10 counts).  
P-values were adjusted using the **Benjamini–Hochberg method**.  
Significant GO terms (padj < 0.05) were **grouped manually into supercategories** using pattern matching (`stringr`).  




