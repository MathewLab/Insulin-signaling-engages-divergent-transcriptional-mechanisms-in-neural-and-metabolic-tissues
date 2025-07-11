# 04_insulin_core_regulation

This module analyzes the differential impact of reduced insulin signaling on core insulin pathway components in OSNs and Fatbody tissues.

## Overview

This figure combines four analytical components:
1. Venn Diagram of significantly regulated insulin core genes (common vs unique).
2. Pairwise Scatter Plot showing directional regulation (log2FC OSNs vs Fatbody).
3. GO Enrichment analysis and Visualisation using `GO_enrichment_circlize.R`.
4. Reactome-style Pathway Schematic highlighting affected genes and regulation direction.

### 1. Venn Diagram – Insulin Core Components (Figure 3A)
This panel uses the same R script as Figure 1 for Venn-Diagram adopted for,a filtered list of insulin core genes.

### 2. Pairwise Scatter Plot – Insulin Core Components (Figure 3B)

This plot shows the direction and strength of log2 fold changes for insulin core genes in OSNs vs Fatbody.
Script modified from figure 1 Scatter Plot Script.

### 3. GO Enrichment – Circos Visualization (Figure 3C)

GO enrichment analysis was performed using **topGO (v2.56.0)** with the **elim algorithm** and **Fisher’s exact test**.  
Each tissue used its own expressed gene background (≥10 counts).  
P-values were adjusted using the **Benjamini–Hochberg method**.  
Significant GO terms (padj < 0.05) were **grouped manually into supercategories** using pattern matching (`stringr`).  
Significant GO terms (padj < 0.05) were **grouped manually into supercategories** using pattern matching (`stringr`). 

### 4. Reactome-style schematic (Figure 3D) 
Canonical insulin signaling pathway genes were visualized using a Reactome-style schematic created with DiagrammeR.  
Node fill colors indicate tissue-specific significance (Red = OSNs, Blue = Fatbody, Purple = Both),  
and directional regulation is shown with triangular glyphs (up/downregulation).




