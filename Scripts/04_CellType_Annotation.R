#===============================================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
#
# Script: 04_CellType_Annotation.R
#
# Purpose:
#   1. Load the Harmony-clustered Seurat object.
#   2. Identify marker genes for each cluster.
#   3. Select the strongest markers for interpretation.
#   4. Check canonical lineage markers.
#   5. Assign broad cell types using DE + canonical markers.
#   6. Visualize and save annotation figures.
#   7. Save the annotated Seurat object.
#===============================================================================


#===============================================================================
# 0. Load packages and clustered object
#===============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

gbm <- readRDS("data/gbm_harmony_clustered.rds")

# Check object
gbm

# Check clusters
table(gbm$seurat_clusters)

# Use Seurat clusters as identities
Idents(gbm) <- "seurat_clusters"

levels(Idents(gbm))

# Check active assay
DefaultAssay(gbm)


#===============================================================================
# 1. Find marker genes for all clusters
#===============================================================================

# Purpose:
# FindAllMarkers() identifies genes that are enriched in each cluster
# relative to the remaining cells.
#
# only.pos = TRUE
#   Keep only genes upregulated in each cluster.
#
# min.pct = 0.25
#   Gene must be detected in at least 25% of cells in either group.
#
# logfc.threshold = 0.25
#   Require at least 0.25 log2 fold-change.

gbm.markers <- FindAllMarkers(
  gbm,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(gbm.markers)


#-------------------------------------------------------------------------------
# Save complete marker table
#-------------------------------------------------------------------------------

write.csv(
  gbm.markers,
  file = "data/gbm_all_cluster_markers.csv",
  row.names = FALSE
)


#===============================================================================
# 2. Select top 10 markers per cluster
#===============================================================================

# Purpose:
# Select the strongest statistically significant markers for easier
# biological interpretation.
#
# avg_log2FC:
#   magnitude of enrichment in the cluster.
#
# pct.1:
#   proportion of cells inside the cluster expressing the gene.
#
# pct.2:
#   proportion of cells outside the cluster expressing the gene.

top10_markers <- gbm.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

top10_markers


#-------------------------------------------------------------------------------
# Save top 10 markers inside project
#-------------------------------------------------------------------------------

write.csv(
  top10_markers,
  file = "data/gbm_top10_markers_per_cluster.csv",
  row.names = FALSE
)


#-------------------------------------------------------------------------------
# Save another copy to Windows Downloads
#-------------------------------------------------------------------------------

write.csv(
  top10_markers,
  file = "C:/Users/Moe9621/Downloads/gbm_top10_markers_per_cluster.csv",
  row.names = FALSE
)


#===============================================================================
# 3. Check canonical markers
#===============================================================================

# Purpose:
#
# Cluster DE:
#   tells us WHAT makes a cluster different.
#
# Canonical markers:
#   help determine WHICH biological lineage the cluster belongs to.
#
# The marker panel covers major cell populations expected in the GBM
# tumor microenvironment.


canonical_markers <- c(
  
  # Myeloid / microglia / macrophages
  "LST1", "TYROBP", "AIF1", "PTPRC",
  
  # T cells
  "CD3D", "CD3E", "TRBC1",
  
  # Mature oligodendrocytes
  "MAG", "MOG", "MOBP", "CLDN11",
  
  # Endothelial
  "PECAM1", "VWF", "CLDN5",
  
  # Fibroblast / stromal
  "COL1A1", "COL1A2", "DCN", "LUM",
  
  # Pericyte / mural
  "RGS5", "PDGFRB", "CSPG4", "MCAM",
  
  # Glial / glioma-associated
  "SOX2", "EGFR", "OLIG2", "PDGFRA",
  
  # Astrocytic program
  "GFAP", "AQP4", "ALDH1L1",
  
  # Cycling
  "MKI67", "TOP2A"
)


#-------------------------------------------------------------------------------
# Canonical marker DotPlot
#-------------------------------------------------------------------------------

# Dot size:
#   percentage of cells in the cluster expressing the gene.
#
# Dot colour:
#   average expression level of the gene in the cluster.

p_dotplot <- DotPlot(
  gbm,
  features = canonical_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Canonical Markers Across GBM Clusters")

p_dotplot


#-------------------------------------------------------------------------------
# Save canonical marker DotPlot
#-------------------------------------------------------------------------------

ggsave(
  filename = "data/gbm_canonical_marker_dotplot.png",
  plot = p_dotplot,
  width = 12,
  height = 8,
  dpi = 300
)


#===============================================================================
# 4. Broad cell-type annotation
#===============================================================================

# Annotation is based on BOTH:
#
# 1. Cluster-specific DE genes from FindAllMarkers()
# 2. Expression of canonical lineage markers
#
# Clusters without convincing normal-cell lineage markers are conservatively
# retained as "Candidate neoplastic".


#-------------------------------------------------------------------------------
# Start with all clusters as candidate neoplastic
#-------------------------------------------------------------------------------

gbm$broad_celltype <- "Candidate neoplastic"


#-------------------------------------------------------------------------------
# Myeloid
#-------------------------------------------------------------------------------
#
# Cluster 2:
#   macrophage/myeloid marker profile.
#
# Cluster 9:
#   P2RY12 and TMEM119 plus strong canonical myeloid markers,
#   supporting microglial identity.
#
# Cluster 14:
#   strong cycling genes such as MKI67, RRM2, TK1 and CEP55,
#   together with LST1, TYROBP, AIF1 and PTPRC.
#
# Therefore cluster 14 represents cycling myeloid cells.

gbm$broad_celltype[
  gbm$seurat_clusters %in% c("2", "9", "14")
] <- "Myeloid"


#-------------------------------------------------------------------------------
# Oligodendrocytes
#-------------------------------------------------------------------------------
#
# Cluster 10:
# MAG, MOG, MOBP, CLDN11 and ERMN
# strongly support mature oligodendrocyte identity.

gbm$broad_celltype[
  gbm$seurat_clusters %in% c("10")
] <- "Oligodendrocytes"


#-------------------------------------------------------------------------------
# T cells
#-------------------------------------------------------------------------------
#
# Cluster 13:
# CD3D, CD3E, CD3G, CD8B and GZMA
# support T-cell identity.

gbm$broad_celltype[
  gbm$seurat_clusters %in% c("13")
] <- "T cells"


#-------------------------------------------------------------------------------
# Fibroblast
#-------------------------------------------------------------------------------
#
# Cluster 15:
# COL1A1, COL3A1, DCN, LUM and CD248
# support a fibroblast/stromal identity.
#
# Some perivascular features are also present, but "Fibroblast" is retained
# as the broad label for comparison with the original paper.

gbm$broad_celltype[
  gbm$seurat_clusters %in% c("15")
] <- "Fibroblast"


#-------------------------------------------------------------------------------
# Endothelial
#-------------------------------------------------------------------------------
#
# Cluster 16:
# VWF, SOX17, ECSCR and canonical vascular markers
# support endothelial identity.

gbm$broad_celltype[
  gbm$seurat_clusters %in% c("16")
] <- "Endothelial"


#===============================================================================
# 5. Check final annotation
#===============================================================================

# Cluster-to-cell-type mapping
table(
  gbm$seurat_clusters,
  gbm$broad_celltype
)

# Number of cells in each broad population
table(gbm$broad_celltype)


#===============================================================================
# 6. Visualize broad cell-type annotation
#===============================================================================

p_annotation <- DimPlot(
  gbm,
  reduction = "umap_after",
  group.by = "broad_celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Broad Cell-Type Annotation")

p_annotation


#-------------------------------------------------------------------------------
# Save annotated UMAP
#-------------------------------------------------------------------------------

ggsave(
  filename = "data/gbm_broad_celltype_annotation.png",
  plot = p_annotation,
  width = 9,
  height = 7,
  dpi = 300
)


#===============================================================================
# 7. Save annotated Seurat object
#===============================================================================

saveRDS(
  gbm,
  file = "data/gbm_broad_celltype_annotated.rds"
)


#===============================================================================
# Result
#===============================================================================
#
# Broad populations identified:
#
#   Candidate neoplastic
#   Myeloid
#   T cells
#   Oligodendrocytes
#   Fibroblast
#   Endothelial
#
# Key saved outputs:
#
# data/gbm_all_cluster_markers.csv
# data/gbm_top10_markers_per_cluster.csv
# data/gbm_canonical_marker_dotplot.png
# data/gbm_broad_celltype_annotation.png
# data/gbm_broad_celltype_annotated.rds
#
# The annotated object can now be used for downstream analyses such as:
#
#   - neoplastic-cell state analysis
#   - myeloid subclustering
#   - cell-type composition analysis
#   - newly diagnosed vs recurrent comparisons
#   - differential expression within defined cell populations
#
#===============================================================================