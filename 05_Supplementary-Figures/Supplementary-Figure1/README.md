# Heatmaps of Top Shared DEGs by Regulatory Pattern

This folder contains heatmaps of the **top 30 shared differentially expressed genes (DEGs)** for each regulatory pattern based on insulin receptor knockdown in *Drosophila melanogaster*:

- **Up–Up**
- **Down–Down**
- **Up–Down**
- **Down–Up**

## Method Summary
Heatmaps were generated using the `pheatmap` R package (v1.0.12).
### Ranking Strategy:
- Genes were ranked within each quadrant by a **combined significance score**:  
  `–log₁₀(padj_OSNs + padj_Fatbody)`
- This scoring method prioritizes genes that are consistently significant across both tissues.

### Heatmap Details:
- **Values visualized**: Log₂ fold change (log₂FC) in each tissue  
- **Clustering**: Row clustering performed using **Euclidean distance**

##  Folder Contents
`scripts/`: R scripts used to generate quadrant-specific heatmaps

  
