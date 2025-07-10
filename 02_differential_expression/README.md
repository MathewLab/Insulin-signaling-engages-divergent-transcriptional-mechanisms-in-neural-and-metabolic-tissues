# Differential Expression Analysis

Differential expression analysis was performed for each dataset using `DESeq2` (v1.44.1) in R (v4.4.1)
by Nevada Bioinformatics Center, University of Nevada, Reno.

## Files
- `DESeq2_results_OSNs_full.csv`: Full differential expression results for OSNs
- `DESeq2_results_FB_full.csv`: Full differential expression results for fat body

DEGs were defined using:
- Adjusted p-value < 0.05 (Benjamini–Hochberg)
- Optional fold-change cutoff: |log₂FC| ≥ 1
