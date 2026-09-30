#===============================================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
# Script: 05_Neoplastic_Cell_Analysis.R
#
# Purpose:
#   This script performs focused analysis of the candidate neoplastic
#   compartment identified during broad cell-type annotation.
#
# Workflow:
#   1. Load the broad-cell-type annotated Seurat object.
#   2. Subset candidate neoplastic cells.
#   3. Re-normalize the neoplastic subset and identify variable features.
#   4. Perform PCA and select informative principal components.
#   5. Assess patient-associated structure before integration.
#   6. Perform Harmony integration across patients.
#   7. Construct the neighborhood graph and identify neoplastic subclusters.
#   8. Visualize the integrated neoplastic-cell landscape.
#   9. Identify cluster-specific marker genes.
#  10. Evaluate published GBM transcriptional-state marker programs.
#  11. Manually annotate AC, MES, NPC and OPC Neftel states.
#  12. Prepare the annotated object for downstream comparison of
#      Initial versus Recurrent tumors.
#
# Input files:
#   - gbm_broad_celltype_annotated.rds
#       Broad-cell-type annotated Seurat object generated in the previous step.
#
#   - GSE229779_cellMetadata.tsv.gz
#       Published cell metadata containing cell IDs, patient information,
#       Initial/Recurrence status, original cell-type annotations,
#       Neftel states and cell-cycle scores.
#
# Output files:
#   - gbm_neoplastic_clustered.rds
#   - gbm_neoplastic_all_markers.csv
#   - gbm_neoplastic_top30_markers_per_cluster.csv
#   - gbm_neoplastic_manual_Neftel_annotated.rds
#
# NOTE:
#   File paths in this script are relative to the project root.
#   Update the paths according to your local project directory structure.
#===============================================================================


#===============================================================================
# 0. Load required packages
#===============================================================================

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)


#===============================================================================
# Load input data
#===============================================================================

# Broad-cell-type annotated Seurat object
gbm <- readRDS(
  "data/gbm_broad_celltype_annotated.rds"
)


# Published cell-level metadata
original_meta <- read.delim(
  "data/GSE229779_cellMetadata (1).tsv.gz",
  check.names = FALSE
)


# Inspect imported metadata
dim(original_meta)
colnames(original_meta)
head(original_meta)
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
# 14. Find markers for every neoplastic subcluster
#
# Differentially expressed markers were identified for each neoplastic
# subcluster using FindAllMarkers().
#
# Initially, the top 10 markers per cluster were inspected. However, for
# several clusters the top 10 genes were dominated by highly specific,
# stress-related, cell-cycle-related, or poorly characterized genes and were
# not sufficient to confidently identify the underlying GBM transcriptional
# state.
#
# Therefore, the top 30 markers were retained to provide a broader view of
# each cluster's transcriptional program and facilitate biological annotation.
#===============================================================================

Idents(gbm_neoplastic) <- "neoplastic_clusters"

neoplastic_markers <- FindAllMarkers(
  gbm_neoplastic,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


#===============================================================================
# 15. Select and save the top 30 markers per cluster
#
# The top 30 significantly enriched genes were ranked by average log2 fold
# change and used as an initial guide for cluster annotation.
#
# Cluster identities were not assigned from individual top markers alone.
# Instead, the top-marker profiles were compared with published marker
# programs defining the major GBM transcriptional states:
#
#   AC  - astrocyte-like
#   MES - mesenchymal-like
#   NPC - neural-progenitor-like
#   OPC - oligodendrocyte-progenitor-like
#
# Coordinated expression of the published marker programs was subsequently
# evaluated using DotPlots. MES1/MES2 and NPC1/NPC2 programs were initially
# examined separately and later collapsed into the four major states for
# downstream analysis.
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


#===============================================================================
# 17. Manual Neftel-state annotation
#
# Neoplastic subclusters were manually annotated by integrating:
#
#   1. The top 30 differentially expressed genes for each cluster
#   2. Published marker programs for AC-, MES-, NPC- and OPC-like GBM states
#   3. DotPlot visualization of coordinated marker expression
#
# MES1/MES2 and NPC1/NPC2 programs were initially evaluated separately,
# then collapsed into the four major Neftel states for downstream analysis.
#
# The original authors' cell-level metadata was additionally inspected as a
# validation reference for selected ambiguous populations.
#
# Clusters without a sufficiently clear AC/MES/NPC/OPC transcriptional
# program were retained as Unassigned.
#===============================================================================

gbm_neoplastic$manual_Neftel_state <- "Unassigned"


# AC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c(
    "0", "4", "6", "12", "13"
  )
] <- "AC"


# MES
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c(
    "1", "9", "11"
  )
] <- "MES"


# NPC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters %in% c(
    "2", "3"
  )
] <- "NPC"


# OPC
gbm_neoplastic$manual_Neftel_state[
  gbm_neoplastic$neoplastic_clusters == "5"
] <- "OPC"


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


# Check available metadata before downstream analysis
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

##===============================================================================
# 22. Differential expression within Neftel states
#
# Differential expression was performed separately within each manually
# annotated Neftel state to compare recurrent versus initial tumors.
#
# ident.1 = Recurrence
# ident.2 = Initial
#
# Therefore:
#   avg_log2FC > 0  = higher expression in Recurrence
#   avg_log2FC < 0  = higher expression in Initial
#
# OPC was not analyzed because no recurrent OPC cells were identified in the
# manual state annotation.
#===============================================================================

ac <- subset(
  gbm_neoplastic,
  subset = manual_Neftel_state == "AC"
)

table(ac$TimePoint)

Idents(ac) <- "TimePoint"

levels(Idents(ac))

DE_AC <- FindMarkers(
  ac,
  ident.1 = "Recurrence",
  ident.2 = "Initial",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

head(DE_AC)

write.csv(
  DE_AC,
  "data/DE_AC_Recurrence_vs_Initial.csv",
  row.names = TRUE
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
  row.names = TRUE
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
  row.names = TRUE
)
#===============================================================================
# 23. Extract significant differentially expressed genes
#
# Differential-expression comparison:
#
#   Recurrence vs Initial
#
# Therefore:
#   avg_log2FC >  0.25 = higher expression in Recurrence
#   avg_log2FC < -0.25 = higher expression in Initial
#
# The terms "higher in Recurrence" and "higher in Initial" are used instead
# of "upregulated" and "downregulated" to make the direction of the comparison
# explicit.
#===============================================================================

#AC
AC_higher_Recurrence <- rownames(
  subset(
    DE_AC,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
)

AC_higher_Initial <- rownames(
  subset(
    DE_AC,
    p_val_adj < 0.05 & avg_log2FC < -0.25
  )
)

# Number of significant genes
length(AC_higher_Recurrence)
length(AC_higher_Initial)

# MES
MES_higher_Recurrence <- rownames(
  subset(
    DE_MES,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
)

MES_higher_Initial <- rownames(
  subset(
    DE_MES,
    p_val_adj < 0.05 & avg_log2FC < -0.25
  )
)

# Number of significant genes
length(MES_higher_Recurrence)
length(MES_higher_Initial)

# NPC
NPC_higher_Recurrence <- rownames(
  subset(
    DE_NPC,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
)

NPC_higher_Initial <- rownames(
  subset(
    DE_NPC,
    p_val_adj < 0.05 & avg_log2FC < -0.25
  )
)

# Number of significant genes
length(NPC_higher_Recurrence)
length(NPC_higher_Initial)

#===============================================================================
# 24. GO Biological Process enrichment
#
# GO enrichment was performed separately for genes with significantly higher
# expression in recurrent and initial tumors within each Neftel state.
#===============================================================================
#================================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

#===============================================================================
# AC state
#===============================================================================

#-------------------------------------------------------------------------------
# AC: genes with higher expression in Recurrence
#-------------------------------------------------------------------------------
GO_AC_Recurrence <- enrichGO(
  gene = AC_higher_Recurrence,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#To check the data
head(as.data.frame(GO_AC_Recurrence))

head(as.data.frame(GO_AC_Recurrence), 15)[,
                                          c("ID", "Description", "GeneRatio", "Count", "p.adjust")
]

# Save AC recurrence-up GO results

write.csv(
  as.data.frame(GO_AC_Recurrence),
  "data/GO_AC_Recurrence.csv"
)

#-------------------------------------------------------------------------------
# AC: genes with higher expression in Initial tumors
#-------------------------------------------------------------------------------

GO_AC_Initial <- enrichGO(
  gene = AC_higher_Initial,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#To check the data
head(as.data.frame(GO_AC_Initial))

write.csv(
  as.data.frame(GO_AC_Initial),
  "data/GO_AC_Initial.csv"
)


#===============================================================================
# MES state
#===============================================================================

# MES: genes with higher expression in Recurrence

GO_MES_Recurrence <- enrichGO(
  gene = MES_higher_Recurrence,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


#to check the data
head(as.data.frame(GO_MES_Recurrence))

write.csv(
  as.data.frame(GO_MES_Recurrence),
  "data/GO_MES_Recurrence.csv"
)

# MES: genes with higher expression in Initial tumors

GO_MES_Initial <- enrichGO(
  gene = MES_higher_Initial,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#to check the data
head(as.data.frame(GO_MES_Initial))

write.csv(
  as.data.frame(GO_MES_Initial),
  "data/GO_MES_down_Initial.csv"
)


#===============================================================================
# NPC state
#===============================================================================

# NPC: genes with higher expression in Recurrence

GO_NPC_Recurrence  <- enrichGO(
  gene = NPC_higher_Recurrence,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#to check the data
head(as.data.frame(GO_NPC_Recurrence))

write.csv(
  as.data.frame(GO_NPC_Recurrence),
  "data/GO_NPC_up_Recurrence.csv",
  row.names = FALSE
)

# NPC: genes with higher expression in Initial tumors
GO_NPC_Initial <- enrichGO(
  gene = NPC_higher_Initial,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#to check the data
head(as.data.frame(GO_NPC_Initial))

write.csv(
  as.data.frame(GO_NPC_Initial),
  "data/GO_NPC_Initial.csv",
  row.names = FALSE
)

#===============================================================================
# Review top GO Biological Processes
#===============================================================================

head(as.data.frame(GO_AC_Initial)[,
                                  c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_AC_Recurrence)[,
                                     c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_MES_Initial)[,
                                   c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_MES_Recurrence)[,
                                      c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_NPC_Initial)[,
                                   c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

head(as.data.frame(GO_NPC_Recurrence)[,
                                      c("ID", "Description", "GeneRatio", "Count", "p.adjust")
], 10)

#===============================================================================
# GO dot plot: AC Recurrence-up genes
#===============================================================================

dotplot(
  GO_AC_Initial,
  showCategory = 15
) +
  ggtitle("AC State - GO BP Enrichment: Higher in Initial")

dotplot(
  GO_AC_Recurrence,
  showCategory = 15
) +
  ggtitle("AC State - GO BP Enrichment: Higher in Recurrence")

dotplot(
  GO_MES_Initial,
  showCategory = 15
) +
  ggtitle("MES State - GO BP Enrichment: Higher in Initial")

dotplot(
  GO_MES_Recurrence,
  showCategory = 15
) +
  ggtitle("MES State - GO BP Enrichment: Higher in Recurrence")

dotplot(
  GO_NPC_Initial,
  showCategory = 15
) +
  ggtitle("NPC State - GO BP Enrichment: Higher in Initial")

dotplot(
  GO_NPC_Recurrence,
  showCategory = 15
) +
  ggtitle("NPC State - GO BP Enrichment: Higher in Recurrence")

#===============================================================================
# GO summary tables
#===============================================================================

GO_AC_Initial_df <- as.data.frame(GO_AC_Initial)
GO_AC_Recurrence_df <- as.data.frame(GO_AC_Recurrence)

GO_MES_Initial_df <- as.data.frame(GO_MES_Initial)
GO_MES_Recurrence_df <- as.data.frame(GO_MES_Recurrence)

GO_NPC_Initial_df <- as.data.frame(GO_NPC_Initial)
GO_NPC_Recurrence_df <- as.data.frame(GO_NPC_Recurrence)

# Save complete GO enrichment tables

write.csv(
  GO_AC_Initial_df,
  "data/GO_AC_Initial.csv",
  row.names = FALSE
)

write.csv(
  GO_AC_Recurrence_df,
  "data/GO_AC_Recurrence.csv",
  row.names = FALSE
)

write.csv(
  GO_MES_Initial_df,
  "data/GO_MES_Initial.csv",
  row.names = FALSE
)

write.csv(
  GO_MES_Recurrence_df,
  "data/GO_MES_Recurrence.csv",
  row.names = FALSE
)

write.csv(
  GO_NPC_Initial_df,
  "data/GO_NPC_Initial.csv",
  row.names = FALSE
)

write.csv(
  GO_NPC_Recurrence_df,
  "data/GO_NPC_Recurrence.csv",
  row.names = FALSE
)

#===============================================================================
# Top 10 GO Biological Processes for each comparison
#
# The 10 most significantly enriched GO Biological Processes were selected
# separately for genes with higher expression in Initial and Recurrence within
# each Neftel state.
#===============================================================================
GO_summary <- bind_rows(
  GO_AC_Initial_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "AC", Comparison = "Higher in Initial"),
  
  GO_AC_Recurrence_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "AC", Comparison = "Higher in Recurrence"),
  
  GO_MES_Initial_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "MES", Comparison = "Higher in Initial"),
  
  GO_MES_Recurrence_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "MES", Comparison = "Higher in Recurrence"),
  
  GO_NPC_Initial_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "NPC", Comparison = "Higher in Initial"),
  
  GO_NPC_Recurrence_df %>%
    filter(p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 10) %>%
    mutate(State = "NPC", Comparison = "Higher in Recurrence")
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
GO_plot_data <- GO_summary

GO_plot_data$log10_padj <- -log10(GO_plot_data$p.adjust)

GO_plot_data$Comparison <- factor(
  paste(GO_plot_data$State, GO_plot_data$Comparison, sep = " - "),
  levels = c(
    "AC - Higher in Recurrence",
    "AC - Higher in Initial",
    "MES - Higher in Recurrence",
    "MES - Higher in Initial",
    "NPC - Higher in Recurrence",
    "NPC - Higher in Initial"
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
# Save combined GO enrichment figure
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
# Save final annotated neoplastic Seurat object
#===============================================================================

saveRDS(
  gbm_neoplastic,
  file = "data/gbm_neoplastic_final.rds"
)
