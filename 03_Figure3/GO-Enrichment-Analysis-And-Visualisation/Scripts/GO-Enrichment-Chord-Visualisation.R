
# Comparative GO Term Enrichment Chord Diagram for OSNs and Fatbody

# Load required libraries
library(tidyverse)
library(circlize
        
# Load and combine data
osns <- read_csv("data/GO_O.csv") %>% mutate(Tissue = "OSNs")
fatbody <- read_csv("data/GO.csv") %>% mutate(Tissue = "Fatbody")

go_combined <- bind_rows(osns, fatbody) %>%
  filter(adj.pval < 0.05)

# Group GO terms by category
group_go_term <- function(term, ontology) {
  term <- tolower(term)
  if (ontology == "MF") {
    if (str_detect(term, "transcription|rna polymerase|dna binding|transcription factor")) "Transcriptional Regulation"
    else if (str_detect(term, "kinase|phosphorylation")) "Kinase Activity"
    else if (str_detect(term, "transporter|channel|pump")) "Transporter Activity"
    else if (str_detect(term, "oxidoreductase|dehydrogenase|redox")) "Oxidoreductase"
    else if (str_detect(term, "atp|gdp|catalytic")) "ATP Binding / Catalysis"
    else "Other (MF)"
  } else if (ontology == "BP") {
    if (str_detect(term, "synapse|neurotransmitter|axon")) "Synaptic Signaling"
    else if (str_detect(term, "metabolic|biosynthetic|catabolic|glycolysis")) "Metabolic Process"
    else if (str_detect(term, "translation|peptide biosynthetic")) "Protein Translation"
    else if (str_detect(term, "cytoskeleton|microtubule|actin")) "Cytoskeleton Dynamics"
    else if (str_detect(term, "vesicle|exocytosis|endocytosis")) "Vesicle Trafficking"
    else if (str_detect(term, "stress|detox|oxidative stress")) "Stress & Detox"
    else "Other (BP)"
  } else if (ontology == "CC") {
    if (str_detect(term, "plasma membrane")) "Plasma Membrane"
    else if (str_detect(term, "cytoskeleton")) "Cytoskeleton"
    else if (str_detect(term, "mitochondrion|respiratory chain")) "Mitochondrion"
    else if (str_detect(term, "vesicle|golgi|endoplasmic reticulum")) "Vesicle & Golgi"
    else if (str_detect(term, "ribosome")) "Ribosome"
    else "Other (CC)"
  } else "Other"
}

# Add supergroups and ontology classes
go_combined <- go_combined %>%
  mutate(
    Supergroup = mapply(group_go_term, Term, Ontology),
    OntoGroup = case_when(
      Supergroup %in% c("Plasma Membrane", "Cytoskeleton", "Mitochondrion", "Vesicle & Golgi", "Ribosome", "Other (CC)") ~ "Cellular Component",
      Supergroup %in% c("Transcriptional Regulation", "Kinase Activity", "Transporter Activity", "Oxidoreductase", "ATP Binding / Catalysis", "Other (MF)") ~ "Molecular Function",
      TRUE ~ "Biological Process"
    )
  )

# Prepare matrix for circos plot
chord_data <- go_combined %>%
  count(Tissue, Supergroup) %>%
  pivot_wider(names_from = Tissue, values_from = n, values_fill = 0) %>%
  column_to_rownames("Supergroup") %>%
  as.matrix() %>%
  t()

# Set sector and color settings

BP_order <- c("Synaptic Signaling", "Protein Translation", "Metabolic Process", "Cytoskeleton Dynamics", "Vesicle Trafficking", "Stress & Detox", "Other (BP)")
CC_order <- c("Plasma Membrane", "Cytoskeleton", "Mitochondrion", "Vesicle & Golgi", "Ribosome", "Other (CC)")
MF_order <- c("Transcriptional Regulation", "Kinase Activity", "Transporter Activity", "Oxidoreductase", "ATP Binding / Catalysis", "Other (MF)")
manual_order <- c("Fatbody", "OSNs", CC_order, MF_order, BP_order)
sector_order <- manual_order[manual_order %in% union(rownames(chord_data), colnames(chord_data))]

group_colors <- c(
  # BP
  "Synaptic Signaling" = "#D73027", "Metabolic Process" = "#FC8D59", "Protein Translation" = "#FEE090",
  "Cytoskeleton Dynamics" = "#cbc9e2", "Vesicle Trafficking" = "#4575B4", "Stress & Detox" = "#E0F3F8", "Other (BP)" = "#999999",
  # MF
  "Transcriptional Regulation" = "#66C2A5", "Kinase Activity" = "#3288BD", "Transporter Activity" = "#ABDDA4",
  "Oxidoreductase" = "#D53E4F", "ATP Binding / Catalysis" = "#FDAE61", "Other (MF)" = "#CCCCCC",
  # CC
  "Plasma Membrane" = "#8DA0CB", "Cytoskeleton" = "#E78AC3", "Mitochondrion" = "#A6D854",
  "Vesicle & Golgi" = "#FFD92F", "Ribosome" = "#E5C494", "Other (CC)" = "#B3B3B3",
  # Tissues
  "Fatbody" = "#a3c8f0", "OSNs" = "#fca9aa"
)

# Chord Diagram Plotting

pdf("results/GO_ChordDiagram_withOuterRings_andLegends.pdf", width = 7, height = 10)

layout(matrix(c(1, 2), nrow = 2), heights = c(7, 3))
par(mar = c(1, 1, 2, 1))  # Circ plot

circos.clear()
circos.par(canvas.xlim = c(-2, 2), canvas.ylim = c(-2, 2), gap.degree = 1, start.degree = 90, track.margin = c(0.01, 0.01))

chordDiagram(
  chord_data,
  grid.col = group_colors,
  order = sector_order,
  transparency = 0.3,
  annotationTrack = "grid",
  annotationTrackHeight = 0.05,
  preAllocateTracks = list(track.height = 0.09)
)

# Outer rings
highlight.sector(CC_order[CC_order %in% sector_order], track.index = 1, col = "#8DA0CB22", border = "#8DA0CB", lwd = 1)
highlight.sector(MF_order[MF_order %in% sector_order], track.index = 1, col = "#66C2A522", border = "#66C2A5", lwd = 1)
highlight.sector(BP_order[BP_order %in% sector_order], track.index = 1, col = "#FC8D5922", border = "#FC8D59", lwd = 1)

circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
  sector <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[2] + 0.13, sector, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
})

title("Comparative GO Term Enrichment in OSNs and Fatbody", cex.main = 1.2, line = -1)

# Add Legends
par(mar = c(0, 0, 0, 0))
plot.new()
legend("topleft",
       legend = c("Fatbody", "OSNs", "Cellular Component", "Molecular Function", "Biological Process"),
       fill = c("#a3c8f0", "#fca9aa", "#8DA0CB", "#66C2A5", "#FC8D59"),
       border = NA, cex = 0.95, title = "Tissue / Ontology", bty = "n")

legend("topright",
       legend = sector_order[!sector_order %in% c("Fatbody", "OSNs")],
       fill = group_colors[sector_order[!sector_order %in% c("Fatbody", "OSNs")]],
       border = NA, cex = 0.85, title = "GO Supergroups", bty = "n")

dev.off()
