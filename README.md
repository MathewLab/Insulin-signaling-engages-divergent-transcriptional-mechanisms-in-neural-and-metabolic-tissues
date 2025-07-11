# Insulin Signaling Drives Tissue-Specific Transcriptomic Programs in Drosophila

This repository contains the full computational workflow, processed data, and figure scripts associated with the manuscript:

**"Insulin signaling engages divergent transcriptional mechanisms in neural and metabolic tissues"**  
*Roshni Jain, Rutuj Kolhe, Cassandra Hui, Juli Petereit, Dennis Mathew*

##  Project Overview

This study presents a comparative transcriptomic analysis of *Drosophila melanogaster* **olfactory sensory neurons (OSNs)** and **fat body (Fb)** under reduced insulin signaling. The findings reveal that insulin signaling coordinates both **shared** and **tissue-specific transcriptional programs**, impacting key pathways related to protein synthesis, synaptic signaling, and metabolism.

##  Data Sources

- **OSNs dataset**: [NCBI BioProject PRJNA1187561](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1187561)  
- **Fat body dataset**: [GEO GSE97447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97447)

## Key Visualizations

- **Figure 1**: DEGs in OSNs and fat body, with shared vs. unique gene overlap and quadrant classification  
- **Figure 2**: KEGG and Reactome pathway enrichment using triangular bubble plots  
- **Figure 3**: Insulin pathway components and GO enrichment (venn + scatter + circos + schematic)  
- **Figure 4**: PPI network of inversely regulated genes (Cytoscape)  
- **Supplementary Figures 1–6**: Regulatory gene families (TFs, kinases, phosphatases, ion channels, GPCRs) across tissues

##  Tools & Packages

- `DESeq2`, `pheatmap`, `ggplot2`, `clusterProfiler`, `ReactomePA`, `circlize`, `VennDiagram`, `DiagrammeR`
- Network analysis using Cytoscape  and stringApp

##  Access and Reproducibility

All code, processed data, and figure scripts are included in this repository. The full analysis is reproducible with R ≥ 4.4.1.

>  [Link to manuscript](#) (to be updated upon publication)

## Contact

For questions or collaboration, please contact:  
**Dr. Dennis Mathew** — dennismathew@unr.edu  
Department of Biology, University of Nevada, Reno

## License

This project is released under the GPL-3.0 License.

