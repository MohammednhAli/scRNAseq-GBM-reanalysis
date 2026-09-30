#===============================================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
#
# Script: 05_Neoplastic_Cell_Analysis.R
#
# Purpose:
#   1. Load the broad-cell-type annotated Seurat object.
#   2. Subset candidate neoplastic cells.
#   3. Recalculate variable features and PCA within the neoplastic compartment.
#   4. Select informative PCs using ElbowPlot and DimHeatmap.
#   5. Visualize patient-associated structure before Harmony.
#   6. Perform Harmony integration on the neoplastic subset.
#   7. Construct the neighborhood graph and neoplastic subclusters.
#   8. Visualize patients and neoplastic subclusters after Harmony.
#   9. Continue with marker identification and neoplastic-state annotation.
#===============================================================================


#===============================================================================
# NEOPLASTIC GBM SUBCLUSTERING PIPELINE
#===============================================================================

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(presto)

gbm <- readRDS("data/gbm_broad_celltype_annotated.rds")

original_meta <- read.delim(
  "data/GSE229779_cellMetadata.tsv.gz"
)
#===============================================================================
# 1. Subset candidate neoplastic cells
#===============================================================================

gbm_neoplastic <- subset(
  gbm,
  subset = broad_celltype == "Candidate neoplastic"
)

gbm_neoplastic


#===============================================================================
# 2. Preserve original broad-cluster identity
#===============================================================================

gbm_neoplastic$parent_cluster <- gbm_neoplastic$seurat_clusters


#===============================================================================
# 3. Normalize
#===============================================================================

gbm_neoplastic <- NormalizeData(
  gbm_neoplastic,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


#===============================================================================
# 4. Identify variable features
#===============================================================================

gbm_neoplastic <- FindVariableFeatures(
  gbm_neoplastic,
  selection.method = "vst",
  nfeatures = 2000
)

#===============================================================================
# 5. Scale data
#===============================================================================

gbm_neoplastic <- ScaleData(
  gbm_neoplastic
)


#===============================================================================
# 6. PCA
#===============================================================================

gbm_neoplastic <- RunPCA(
  gbm_neoplastic,
  features = VariableFeatures(gbm_neoplastic),
  reduction.name = "pca_neoplastic"
)



#===============================================================================
# 7. PC selection
#===============================================================================

ElbowPlot(
  gbm_neoplastic,
  reduction = "pca_neoplastic"
)

DimHeatmap(
  gbm_neoplastic,
  dims = 1:10,
  cells = 500,
  balanced = TRUE,
  reduction = "pca_neoplastic"
)

# PCs 1:10 selected


#===============================================================================
# 8. UMAP before integration
#===============================================================================

gbm_neoplastic <- RunUMAP(
  gbm_neoplastic,
  reduction = "pca_neoplastic",
  dims = 1:10,
  reduction.name = "umap_neoplastic_before"
)

DimPlot(
  gbm_neoplastic,
  reduction = "umap_neoplastic_before",
  group.by = "orig.ident"
) +
  ggtitle("Neoplastic Cells Before Harmony")


#===============================================================================
# 9. Harmony integration
#===============================================================================

gbm_neoplastic <- RunHarmony(
  gbm_neoplastic,
  "orig.ident",
  reduction.use = "pca_neoplastic",
  reduction.save = "harmony_neoplastic"
)

#===============================================================================
# 10. Find neighbors
#===============================================================================

gbm_neoplastic <- FindNeighbors(
  gbm_neoplastic,
  reduction = "harmony_neoplastic",
  dims = 1:10
)

#===============================================================================
# 11. Cluster neoplastic cells
#===============================================================================

gbm_neoplastic <- FindClusters(
  gbm_neoplastic,
  resolution = 0.5
)

gbm_neoplastic$neoplastic_clusters <- gbm_neoplastic$seurat_clusters

table(gbm_neoplastic$neoplastic_clusters)
#===============================================================================
# 12. UMAP after Harmony
#===============================================================================

gbm_neoplastic <- RunUMAP(
  gbm_neoplastic,
  reduction = "harmony_neoplastic",
  dims = 1:10,
  reduction.name = "umap_neoplastic_after"
)

DimPlot(
  gbm_neoplastic,
  reduction = "umap_neoplastic_after",
  group.by = "neoplastic_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Neoplastic Subclusters")

#===============================================================================
# 13. Save clustered object
#===============================================================================

saveRDS(
  gbm_neoplastic,
  file = "data/gbm_neoplastic_clustered.rds"
)

#===============================================================================
# 14. Find markers for every neoplastic cluster
#===============================================================================

Idents(gbm_neoplastic) <- "neoplastic_clusters"

neoplastic_markers <- FindAllMarkers(
  gbm_neoplastic,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#===============================================================================
# 14. Find markers for every neoplastic cluster
#===============================================================================

Idents(gbm_neoplastic) <- "neoplastic_clusters"

neoplastic_markers <- FindAllMarkers(
  gbm_neoplastic,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#===============================================================================
# 16. Save marker tables
#===============================================================================

write.csv(
  neoplastic_markers,
  "data/gbm_neoplastic_all_markers.csv",
  row.names = FALSE
)

top30_neoplastic <- neoplastic_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 30,
    with_ties = FALSE
  ) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  top30_neoplastic,
  "data/gbm_neoplastic_top30_markers_per_cluster.csv",
  row.names = FALSE
)




##===============================================================================
# 17. Manual Neftel-state annotation
#
# Neoplastic clusters were annotated using:
#   1. Top differentially expressed genes for each cluster
#   2. Published marker programs for AC-, MES-, NPC- and OPC-like GBM states
#   3. DotPlot visualization of coordinated marker expression
#
# MES1/MES2 and NPC1/NPC2 programs were initially evaluated separately,
# then collapsed into the four major Neftel states for downstream analysis.
#
# Manual annotations were subsequently compared with the original authors'
# cell-level Neftel_State metadata as validation.


gbm_neoplastic$manual_Neftel_state <- "Unassigned"

# AC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c("0", "4", "6", "12", "13")
] <- "AC"

# MES
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c("1", "9", "11")
] <- "MES"

# NPC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c("2", "3")
] <- "NPC"

# OPC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters == "5"
] <- "OPC"

# NPC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c(
    "2", "3"
  )
] <- "NPC"


#===============================================================================
# 18. Check manual annotations
#===============================================================================

table(
  gbm_neoplastic$neoplastic_clusters,
  gbm_neoplastic$manual_Neftel_state
)


#===============================================================================
# 19. Visualize neoplastic subclusters and manual Neftel states
#===============================================================================

p_clusters <- DimPlot(
  gbm_neoplastic,
  reduction = "umap_neoplastic_after",
  group.by = "neoplastic_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Neoplastic Subclusters")


p_states <- DimPlot(
  gbm_neoplastic,
  reduction = "umap_neoplastic_after",
  group.by = "manual_Neftel_state",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Neftel-State Annotation")


p_neoplastic_annotation <- p_clusters + p_states

p_neoplastic_annotation

colnames(gbm_neoplastic@meta.data)



#===============================================================================
# 20. Add Initial / Recurrence information
#===============================================================================

idx <- match(
  colnames(gbm_neoplastic),
  original_meta$cellID
)

# Check all cells matched
sum(is.na(idx))

gbm_neoplastic$TimePoint <-
  original_meta$mergeTimePoint[idx]


table(gbm_neoplastic$TimePoint)

table(
  gbm_neoplastic$orig.ident,
  gbm_neoplastic$TimePoint)

#===============================================================================
# 21. Save final neoplastic annotated object
#===============================================================================

saveRDS(
  gbm_neoplastic,
  file = "data/gbm_neoplastic_annotated.rds"
)

# 22. Differential expression: AC state
#===============================================================================

ac <- subset(
  gbm_neoplastic,
  subset = manual_Neftel_state == "AC"
)

table(ac$TimePoint)

Idents(ac) <- "TimePoint"

DE_AC <- FindMarkers(
  ac,
  ident.1 = "Recurrence",
  ident.2 = "Initial",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

write.csv(
  DE_AC,
  "data/DE_AC_Recurrence_vs_Initial.csv",
  row.names = FALSE
)

# Differential expression: MES state
#===============================================================================

mes <- subset(
  gbm_neoplastic,
  subset = manual_Neftel_state == "MES"
)

table(mes$TimePoint)

Idents(mes) <- "TimePoint"

DE_MES <- FindMarkers(
  mes,
  ident.1 = "Recurrence",
  ident.2 = "Initial",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

write.csv(
  DE_MES,
  "data/DE_MES_Recurrence_vs_Initial.csv",
  row.names = FALSE
)
# Differential expression: NPC state
#===============================================================================

npc <- subset(
  gbm_neoplastic,
  subset = manual_Neftel_state == "NPC"
)

table(npc$TimePoint)

Idents(npc) <- "TimePoint"

DE_NPC <- FindMarkers(
  npc,
  ident.1 = "Recurrence",
  ident.2 = "Initial",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

write.csv(
  DE_NPC,
  "data/DE_NPC_Recurrence_vs_Initial.csv",
  row.names = FALSE
)
# 23. Extract significant DEGs
#===============================================================================

AC_up <- rownames(
  subset(
    DE_AC,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
)

AC_down <- rownames(
  subset(
    DE_AC,
    p_val_adj < 0.05 & avg_log2FC < -0.25
  )
)

length(AC_up)
length(AC_down)

MES_up <- rownames(
  subset(
    DE_MES,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
)

MES_down <- rownames(
  subset(
    DE_MES,
    p_val_adj < 0.05 & avg_log2FC < -0.25
  )
)

length(MES_up)
length(MES_down)

NPC_up <- rownames(
  subset(
    DE_NPC,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
)

NPC_down <- rownames(
  subset(
    DE_NPC,
    p_val_adj < 0.05 & avg_log2FC < -0.25
  )
)

length(NPC_up)
length(NPC_down)

#===============================================================================
# GO Biological Process Enrichment
#===============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

GO_AC_up <- enrichGO(
  gene = AC_up,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(GO_AC_up))

head(as.data.frame(GO_AC_up), 15)[,
                                  c("ID", "Description", "GeneRatio", "Count", "p.adjust")
]

# Save AC recurrence-up GO results

write.csv(
  as.data.frame(GO_AC_up),
  "data/GO_AC_up_Recurrence.csv",
  row.names = FALSE
)

# AC: Initial-up genes
GO_AC_down <- enrichGO(
  gene = AC_down,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(GO_AC_down),
  "data/GO_AC_down_Initial.csv",
  row.names = FALSE
)

# MES: Recurrence-up genes
GO_MES_up <- enrichGO(
  gene = MES_up,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(GO_MES_up),
  "data/GO_MES_up_Recurrence.csv",
  row.names = FALSE
)

# MES: Initial-up genes
GO_MES_down <- enrichGO(
  gene = MES_down,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(GO_MES_down),
  "data/GO_MES_down_Initial.csv",
  row.names = FALSE
)

# NPC: Recurrence-up genes
GO_NPC_up <- enrichGO(
  gene = NPC_up,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(GO_NPC_up),
  "data/GO_NPC_up_Recurrence.csv",
  row.names = FALSE
)

# NPC: Initial-up genes
GO_NPC_down <- enrichGO(
  gene = NPC_down,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(GO_NPC_down),
  "data/GO_NPC_down_Initial.csv",
  row.names = FALSE
)

#===============================================================================
# Review top GO Biological Processes
#===============================================================================

head(as.data.frame(GO_AC_up)[,
                             c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_AC_down)[,
                               c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_MES_up)[,
                              c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_MES_down)[,
                                c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_NPC_up)[,
                              c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_NPC_down)[,
                                c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

#===============================================================================
# GO dot plot: AC Recurrence-up genes
#===============================================================================

dotplot(
  GO_AC_up,
  showCategory = 15
) +
  ggtitle("AC State - Recurrence-up GO Biological Processes")

dotplot(
  GO_AC_down,
  showCategory = 15
) +
  ggtitle("AC State - Initial-up GO Biological Processes")

dotplot(
  GO_MES_up,
  showCategory = 15
) +
  ggtitle("MES State - Recurrence-up GO Biological Processes")

dotplot(
  GO_MES_down,
  showCategory = 15
) +
  ggtitle("MES State - Initial-up GO Biological Processes")

dotplot(
  GO_NPC_up,
  showCategory = 15
) +
  ggtitle("NPC State - Recurrence-up GO Biological Processes")

dotplot(
  GO_NPC_down,
  showCategory = 15
) +
  ggtitle("NPC State - Initial-up GO Biological Processes")

#===============================================================================
# GO summary tables
#===============================================================================

GO_AC_up_df <- as.data.frame(GO_AC_up)
GO_AC_down_df <- as.data.frame(GO_AC_down)
GO_MES_up_df <- as.data.frame(GO_MES_up)
GO_MES_down_df <- as.data.frame(GO_MES_down)
GO_NPC_up_df <- as.data.frame(GO_NPC_up)
GO_NPC_down_df <- as.data.frame(GO_NPC_down)

write.csv(
  GO_AC_up_df,
  "data/GO_AC_up_Recurrence.csv",
  row.names = FALSE
)

write.csv(
  GO_AC_down_df,
  "data/GO_AC_down_Initial.csv",
  row.names = FALSE
)

write.csv(
  GO_MES_up_df,
  "data/GO_MES_up_Recurrence.csv",
  row.names = FALSE
)

write.csv(
  GO_MES_down_df,
  "data/GO_MES_down_Initial.csv",
  row.names = FALSE
)

write.csv(
  GO_NPC_up_df,
  "data/GO_NPC_up_Recurrence.csv",
  row.names = FALSE
)

write.csv(
  GO_NPC_down_df,
  "data/GO_NPC_down_Initial.csv",
  row.names = FALSE
)

#===============================================================================
# Top 10 GO Biological Processes for each comparison
#===============================================================================

GO_summary <- bind_rows(
  GO_AC_up_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "AC", Comparison = "Recurrence-up"),
  
  GO_AC_down_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "AC", Comparison = "Initial-up"),
  
  GO_MES_up_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "MES", Comparison = "Recurrence-up"),
  
  GO_MES_down_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "MES", Comparison = "Initial-up"),
  
  GO_NPC_up_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "NPC", Comparison = "Recurrence-up"),
  
  GO_NPC_down_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "NPC", Comparison = "Initial-up")
)

write.csv(
  GO_summary,
  "data/GO_summary_top10.csv",
  row.names = FALSE
)

GO_summary[, c(
  "State",
  "Comparison",
  "Description",
  "GeneRatio",
  "Count",
  "p.adjust"
)]

#===============================================================================
# Combined GO enrichment figure
#===============================================================================

GO_plot_data <- GO_summary

GO_plot_data$log10_padj <- -log10(GO_plot_data$p.adjust)

GO_plot_data$Comparison <- factor(
  paste(GO_plot_data$State, GO_plot_data$Comparison, sep = " - "),
  levels = c(
    "AC - Recurrence-up",
    "AC - Initial-up",
    "MES - Recurrence-up",
    "MES - Initial-up",
    "NPC - Recurrence-up",
    "NPC - Initial-up"
  )
)

p_GO_summary <- ggplot(
  GO_plot_data,
  aes(
    x = Comparison,
    y = reorder(Description, log10_padj),
    size = Count,
    color = log10_padj
  )
) +
  geom_point() +
  theme_bw() +
  labs(
    title = "GO Biological Process Enrichment",
    x = "State and comparison",
    y = "Biological Process",
    size = "Gene count",
    color = "-log10 adjusted p-value"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 8
    )
  )

p_GO_summary

#===============================================================================
# Save combined GO figure
#===============================================================================

ggsave(
  "data/GO_summary_combined.pdf",
  plot = p_GO_summary,
  width = 12,
  height = 8
)

ggsave(
  "data/GO_summary_combined.png",
  plot = p_GO_summary,
  width = 12,
  height = 8,
  dpi = 300
)

#===============================================================================
# Save final neoplastic Seurat object
#===============================================================================

saveRDS(
  gbm_neoplastic,
  file = "data/gbm_neoplastic_final.rds"
)
