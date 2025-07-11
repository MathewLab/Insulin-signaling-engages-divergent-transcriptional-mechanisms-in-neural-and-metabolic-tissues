# Supplementary Figures 2–6 – Tissue-Specific Deployment of Gene Families

This folder contains the analysis and visualizations for Supplementary Figures 2 through 6, highlighting tissue-specific regulation of transcription factors, kinases, phosphatases, ion channels, and GPCRs in *Drosophila melanogaster* following insulin receptor knockdow

##  Overview of Figures

Each figure includes three components:
- **(A) Venn Diagram** – Identifies genes differentially expressed in OSNs, fat body, or both (adjusted p < 0.05).
- **(B) Scatter Plot** – Visualizes log₂ fold change of shared genes across tissues.
- **(C) Circos Plot** – Depicts gene family representation and tissue specificity using a chord diagram.

Colors indicate quadrant-based regulation patterns in scatter plots:
- **Red**: Up in OSNs and Fat body (Up–Up)  
- **Blue**: Down in both tissues (Down–Down)  
- **Orange**: Up in OSNs, Down in Fat body (Up–Down)  
- **Purple**: Down in OSNs, Up in Fat body (Down–Up)


##  Figures and Corresponding Scripts

###  Supplementary Figure 2 – Transcription Factors
- **Venn Diagram**: `01_Figure1/Differential-Expression_Analysis_and_Visualisation/Scripts/Venn-Diagram.R`
- **Scatter Plot**: `01_Figure1/Differential-Expression_Analysis_and_Visualisation/Scripts/Scatter-Plot.R`
- **Circos Plot**: `03_Figure3/GO-Enrichment-Analysis-And-Visualisation/Scripts/GO-Enrichment-Chord-Visualisation.R`

###  Supplementary Figure 3 – Kinases  
Same scripts as above; data subset for kinases.

###  Supplementary Figure 4 – Phosphatases  
Same scripts as above; data subset for phosphatases.

###  Supplementary Figure 5 – Ion Channels  
Same scripts as above; data subset for ion channels.

###  Supplementary Figure 6 – GPCRs  
Same scripts as above; data subset for GPCRs.

##  Notes

- Venn and scatter plots were generated using shared scripts across all gene families with customized input tables.
- Circos plots use gene family annotations to visualize tissue-specific deployment with ribbon overlaps reflecting shared subfamilies.




