# plot_enrichment_bubbles.R

library(ggplot2)
library(readr)
library(dplyr)

# Load final filtered file
df <- read_csv("results/enrichment_tables/Comprehensive_Grouped_Pathways.csv")

# Clean pathway labels
df$PathwayShort <- gsub("KEGG_|Reactome_|Pathway_|Drosophila_", "", df$Description)
df$PathwayShort <- substr(df$PathwayShort, 1, 40)
df$PathwayLabel <- df$PathwayShort

# Assign triangle shapes
df$Shape <- ifelse(df$Regulation == "Up", 24, 25)

# Define color palette by group
group_colors <- c(
  "Cell Structure / Adhesion / ECM" = "#E1AB00",
  "Immunity / Inflammation" = "#297373",
  "Metabolic Signaling" = "#F6D5B6",
  "Metabolism: Lipids & Miscellaneous" = "#1696d2",
  "Mitochondrial Function & Energy Metabolism" = "#7FC97F",
  "Neurodegeneration" = "#b5a6d7",
  "Other / Miscellaneous" = "#cfcfcf",
  "Protein Degradation / Turnover" = "#FFA69E",
  "Protein Synthesis / Translation" = "#C4D8A7",
  "Synaptic Signaling" = "#e0783F"
)

# Set facet order (OSNs left, Fatbody right)
df$Dataset <- factor(df$Dataset, levels = c("OSNs", "Fatbody"))

# Select top 20 pathways per group
df_top <- df %>%
  group_by(Dataset, Regulation) %>%
  slice_min(order_by = p.adjust, n = 20, with_ties = FALSE) %>%
  ungroup()

# Order y-axis by significance
df_top$PathwayLabel <- factor(
  df_top$PathwayLabel,
  levels = unique(df_top[order(df_top$Dataset, df_top$Regulation, df_top$p.adjust), "PathwayLabel"][[1]])
)

# Plot
p <- ggplot(df_top, aes(x = Regulation, y = PathwayLabel)) +
  geom_point(
    aes(size = -log10(p.adjust), fill = Group, shape = as.factor(Shape)),
    stroke = 0.2, alpha = 0.93
  ) +
  geom_vline(xintercept = 1:3, linetype = "solid", color = "#F0F0F0", size = 0.4) +
  scale_shape_manual(
    name = "Regulation",
    values = c("24" = 24, "25" = 25),
    labels = c("Upregulated", "Downregulated")
  ) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_size(name = "-log10(p.adjust)", range = c(3, 8)) +
  facet_grid(. ~ Dataset, scales = "free_y", space = "free") +
  labs(
    title = "Convergent Pathway Signatures in Neuronal and Metabolic Tissues Under Reduced Insulin Signaling",
    x = "Regulation",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 9, hjust = 1),
    axis.text.x = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "#F0F0F0", size = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    size = guide_legend(order = 3),
    shape = guide_legend(order = 2)
  )

# Save
ggsave("results/enrichment_tables/Core-Insulin_Signaling_Enrichment_Top20_TRIANGLE_BUBBLE.pdf",
       p, width = 15, height = 11, device = cairo_pdf)

# Show plot
print(p)
