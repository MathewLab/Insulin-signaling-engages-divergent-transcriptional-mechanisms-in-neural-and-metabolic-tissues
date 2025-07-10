# Figure 1 – Transcriptome Comparison Overview

This folder contains the complete analysis and visualizations for Figure 1 of the study,
which compares gene expression between OSNs and fat body following insulin receptor knockdown.

## Panels

- Figure 1A–1B: Volcano plots of DEGs in OSNs and Fatbody
- Figure 1C: Venn diagrams showing gene overlaps
- Figure 1D: Pairwise expression correlation of common DEGs (padj < 0.05)

## Folder Structure

- `scripts/`: R scripts for each sub-panel
- `data/`: Processed tables (e.g., DEG overlap, quadrant classification)

## Notes

- The volcano plots use a shared script (`volcano_plot.R`) with different input files for OSNs vs. fat body.
- Venn diagrams are generated for different p-value and fold-change thresholds.
- Pairwise plots include Pearson correlation and classification into Up-Up, Down-Down, Up-Down, Down-Up quadrants.
