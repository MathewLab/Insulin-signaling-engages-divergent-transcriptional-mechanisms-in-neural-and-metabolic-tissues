# Preprocessing and Quality Control

This folder contains documentation for the preprocessing and QC of RNA-seq data performed by the Nevada Bioinformatics Center.

Scripts used in this pipeline was developed by Nevada Bioinformatics Center (NBC), at the University of Nevada, Reno  and  are **not publicly available** but can be provided **upon request**.

## Steps Performed

1. Adapter trimming with `fastp` (v0.20.0)
2. Quality assessment using `FastQC` and `MultiQC`
3. Alignment with `STAR` (v2.7) to FlyBase r6.60 genome
4. Gene quantification using `featureCounts` (v2.0.0)

## Summary File

- `all.star.fastp_se.fixcol2.xlsx`: Summary of alignment and quality metrics for all RNA-seq samples (provided by NBC).


