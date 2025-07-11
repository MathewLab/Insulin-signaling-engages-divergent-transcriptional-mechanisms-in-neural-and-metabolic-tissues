# Figure 1 – Transcriptome Comparison Overview

This folder contains the complete analysis and visualizations for **Figure 1** of the study, which compares gene expression between **olfactory sensory neurons (OSNs)** and **fat body** following insulin receptor knockdown in *Drosophila melanogaster*.


##  Differential Expression Analysis

Differential expression analysis was performed using `DESeq2` (v1.44.1) in R (v4.4.1) by the **Nevada Bioinformatics Center**, University of Nevada, Reno.

### DEGs were defined using:
- Adjusted p-value < 0.05 (Benjamini–Hochberg correction)
- Optional log₂ fold-change threshold: |log₂FC| ≥ 1

### Input/Output Files
- `InR_OSNs_All .csv – Full results for OSNs
- `InR_Fatbody_All .csv – Full results for fat body

##  Figure 1 Panels

- **Figure 1A & 1B**: Volcano plots of DEGs in OSNs and Fat Body
- **Figure 1C**: Venn diagrams showing gene overlaps across tissues
- **Figure 1D**: Pairwise correlation plot of commonly significant genes (padj < 0.05), categorized by quadrant:
  - Up-Up, Down-Down, Up-Down, Down-Up

## Folder Structure
All analyses and visualizations for Figure 1 are organized under:
`Differential-Expression_Analysis_and_Visualisation/DEG_Results/`
- `scripts/`  
  Contains R scripts used to generate each figure panel (volcano plots, Venn diagrams, pairwise scatter plots).
- `data/`  
  Includes processed result table such as: 
  - Quadrant classification of common DEGs

##  Notes

- The volcano plots are generated using a shared script (`volcano_plot.R`) with distinct inputs for each tissue.
- Venn diagrams are calculated using combinations of adjusted p-value and optional fold-change cutoffs.
- Scatterplots classify common DEGs into quadrants to highlight regulation directionality.



