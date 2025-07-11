# Raw Data Overview, Preprocessing, and Quality Control

This document summarizes the RNA-seq datasets and the preprocessing and quality control steps used in the analysis of olfactory sensory neurons (OSNs) and fat body tissues from *Drosophila melanogaster*.

## BioProject Accessions

- **OSNs Dataset**: [NCBI BioProject PRJNA1187561](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1187561)
- **Fat Body Dataset**: [GEO Accession GSE97447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97447)

These datasets were used to investigate tissue-specific transcriptomic responses under altered insulin signaling conditions.

## Preprocessing and Quality Control

All preprocessing and quality control steps were performed by the Nevada Bioinformatics Center (NBC) at the University of Nevada, Reno.

> **Note:** Scripts used in this pipeline were developed by NBC.

### Steps Performed

1. Adapter trimming using `fastp` (v0.20.0)
2. Quality assessment using `FastQC` and `MultiQC`
3. Alignment using `STAR` (v2.7) to the FlyBase r6.60 genome
4. Gene quantification using `featureCounts` (v2.0.0)

## Summary File

- `all.star.fastp_se.fixcol2.xlsx`: Summary of alignment and quality metrics for all RNA-seq samples (provided by NBC).

