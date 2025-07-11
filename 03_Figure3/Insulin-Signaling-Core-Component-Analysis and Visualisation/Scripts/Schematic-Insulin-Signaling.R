# Description: Generates a Reactome-style schematic for insulin signaling genes,
# showing significance and direction of regulation across OSNs and Fatbody.

# Load libraries
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(readr)
library(pdftools)

# ---- Load and prepare input ----
core <- read_csv("../data/insulin_pathway_gene_list.csv")

# Clean regulation annotations
for (col in c("Regulation OSNs", "Regulation Fatbody")) {
  core[[col]] <- trimws(core[[col]])
  core[[col]][tolower(core[[col]]) == "up"] <- "Up"
  core[[col]][tolower(core[[col]]) == "down"] <- "Down"
  core[[col]][tolower(core[[col]]) == "none"] <- "None"
}

# Assign fill color by tissue-specific significance
core$group <- with(core, ifelse(padj_OSNs < 0.05 & padj_Fatbody < 0.05, "Both",
                                ifelse(padj_OSNs < 0.05, "OSNs",
                                       ifelse(padj_Fatbody < 0.05, "Fatbody", "None"))))
core$color <- ifelse(core$group == "Both", "#bda7e7",
                     ifelse(core$group == "OSNs", "#fca9aa",
                            ifelse(core$group == "Fatbody", "#a3c8f0", "#e6e6e6")))

# Border color = combined directional regulation
core$border_color <- with(core,
  ifelse(`Regulation OSNs` == "Up" & `Regulation Fatbody` == "Up", "#F33A6A",
  ifelse(`Regulation OSNs` == "Down" & `Regulation Fatbody` == "Down", "#2979FF",
  ifelse(`Regulation OSNs` == "Up" & `Regulation Fatbody` == "Down", "#E6B800",
  ifelse(`Regulation OSNs` == "Down" & `Regulation Fatbody` == "Up", "#8BC34A", "gray70")))))

# Create mapping
color_map   <- setNames(core$color, core$gene_symbol)
border_map  <- setNames(core$border_color, core$gene_symbol)

# Core genes in order
core_genes <- c("Ilp2", "Ilp3", "Ilp5", "Ilp6", "Ilp7", "InR", "chico", "Pi3K92E",
                "Pdk1", "Akt", "foxo", "Tsc1", "Rheb", "mTor", "raptor",
                "S6k", "Thor", "sgg", "CycD", "Myc", "eIF4E")
core_genes <- core_genes[core_genes %in% core$gene_symbol]

# Triangular glyphs for regulation direction
osn_tri <- function(gene) {
  reg <- core[core$gene_symbol == gene, ]$`Regulation OSNs`
  if (reg == "Up") return('<font color="#E84A5F">&#9650;</font>')
  if (reg == "Down") return('<font color="#2274A5">&#9660;</font>')
  return("")
}
fb_tri <- function(gene) {
  reg <- core[core$gene_symbol == gene, ]$`Regulation Fatbody`
  if (reg == "Up") return('<font color="#FF914D">&#9650;</font>')
  if (reg == "Down") return('<font color="#009688">&#9660;</font>')
  return("")
}

# Generate DiagrammeR node strings
node_fmt <- function(gene) {
  clr <- color_map[gene]
  bdr <- border_map[gene]
  gene_label <- if (gene == "sgg") "sgg(GSK3-β)" else gene
  label <- sprintf('<<TABLE BORDER="0" CELLBORDER="0">
    <TR><TD><B><FONT FACE="Arial">%s</FONT></B></TD></TR>
    <TR><TD><FONT FACE="Arial">%s&nbsp;%s</FONT></TD></TR>
  </TABLE>>', gene_label, osn_tri(gene), fb_tri(gene))
  sprintf('%s [label = %s, fillcolor="%s", color="%s", fontname="Arial", fontcolor=black,
          shape=box, style="rounded,filled", width=2.0, height=0.7, fixedsize=true, penwidth=3]',
          gene, label, clr, bdr)
}
nodes_code <- paste(sapply(core_genes, node_fmt), collapse = "\n  ")

# Annotated labels
annotations <- '
  neuro [label = "Neurodegeneration", shape=plaintext, fontcolor="#7B4397", fontsize=16]
  protein [label = "Protein Synthesis", shape=plaintext, fontcolor="#5E865E", fontsize=16]
  metabneuro [label = "Metabolism / Neurodegeneration", shape=plaintext, fontcolor="#8b62b9", fontsize=16]
  signaling [label = "Metabolic Signaling", shape=plaintext, fontcolor="#C17D30", fontsize=16]
  immunity [label = "Immunity / Inflammation", shape=plaintext, fontcolor="#297373", fontsize=16]
'

# Edges (pathway structure)
edges_code <- '
  Ilp2 -> InR
  Ilp3 -> InR
  Ilp5 -> InR
  Ilp6 -> InR
  Ilp7 -> InR
  InR -> chico
  chico -> Pi3K92E
  Pi3K92E -> Pdk1
  Pdk1 -> Akt
  Akt -> foxo [label="inhibits", style=dashed, color=gray40]
  Akt -> Tsc1 [label="inhibits", style=dashed, color=gray40]
  Tsc1 -> Rheb [arrowhead=tee, color=gray40]
  Rheb -> mTor
  mTor -> raptor
  raptor -> S6k
  raptor -> Thor
  raptor -> sgg
  raptor -> CycD
  S6k -> Myc
  Thor -> eIF4E

  neuro -> mTor [color="#b5a6d7", style=dashed, arrowhead=none]
  protein -> S6k
  protein -> Thor
  protein -> Myc
  protein -> eIF4E
  metabneuro -> sgg
  metabneuro -> foxo
  signaling -> Akt
  signaling -> Pi3K92E
  immunity -> foxo
'
g <- grViz(sprintf('
digraph insulin_signaling {
  graph [rankdir = TB, bgcolor=white]
  node [fontsize=16, fontname="Arial"]
  %s
  %s
  %s
}', nodes_code, annotations, edges_code))

# Save to file
svg <- export_svg(g)
rsvg_pdf(charToRaw(svg), file = "../plots/page1_pathway.pdf")
writeLines(svg, "../plots/insulin_signaling_pathway.svg")

pdf("../plots/page2_legend.pdf", width=12, height=7)
par(mfrow=c(3,1), mar=c(0,0,0,0))
plot.new()
legend("center", legend = c("OSNs only", "Fatbody only", "Both tissues", "Not significant"),
       fill = c("#fca9aa", "#a3c8f0", "#bda7e7", "#e6e6e6"), border = "gray60",
       title = "Node Fill: Significance by Tissue", horiz = TRUE, cex=1.4, bty="n")
plot.new()
legend("center", legend = c("Up/Up", "Down/Down", "Up(OSNs)/Down(FB)", "Down(OSNs)/Up(FB)", "Mixed/NS"),
       lwd=6, col = c("#F33A6A", "#2979FF", "#E6B800", "#8BC34A", "gray70"),
       title = "Node Border: Combined Regulation", horiz = TRUE, cex=1.4, bty="n")
plot.new()
legend("center", legend = c("OSNs Up", "OSNs Down", "Fatbody Up", "Fatbody Down"),
       pch = c(24,25,24,25), pt.bg=c("#E84A5F", "#2274A5", "#FF914D", "#009688"),
       col=c("#E84A5F", "#2274A5", "#FF914D", "#009688"), pt.cex=2.1,
       title = "Triangles: Regulation by Tissue", horiz=TRUE, cex=1.25, bty="n")
dev.off()

pdf_combine(c("../plots/page1_pathway.pdf", "../plots/page2_legend.pdf"),
            output = "../plots/Reactome_Style_Pathway.pdf")
