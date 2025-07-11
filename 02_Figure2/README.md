# Figure 2: Pathway Enrichment Analysis

This folder contains scripts and results for pathway enrichment analysis of differentially expressed genes in OSNs and Fatbody under reduced insulin signaling.

##  Scripts

Located in `scripts/`:

- `KEGG_enrichment.R`: KEGG-based enrichment using `clusterProfiler` and `msigdbr`.
- `Reactome_enrichment.R`: Reactome enrichment using `ReactomePA`.
- `MSigDB_enrichment.R`: Enrichment using curated gene sets from MSigDB.

Each script:
- Separates analysis for upregulated and downregulated genes.
- Filters significant terms (padj < 0.05).
- Adds mapped gene symbols to output.
- Removes redundant pathway terms using Jaccard-like string similarity.

##  Results
Stored in `results/enrichment_tables/`:
- `Comprehensive_Grouped_Pathways.csv`: Final table used for plotting grouped pathways in Figure 2.

##  Plot Summary

A dot plot was generated showing the top 20 enriched pathways for each tissue and regulation pattern. Pathways are grouped functionally (e.g., Metabolism, Synaptic Signaling), with triangle shapes indicating direction (up/downregulated).

See the plotting script inside `04_figure2_pathway_enrichment/scripts/` for full details.


## Related Tools
- [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
- [`ReactomePA`](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html)
- [`msigdbr`](https://cran.r-project.org/web/packages/msigdbr/)

