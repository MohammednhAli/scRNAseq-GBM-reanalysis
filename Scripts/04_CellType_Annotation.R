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
library(presto)

# 'presto' provides a faster implementation of the Wilcoxon Rank Sum Test
# used by Seurat for FindMarkers() and FindAllMarkers().
# Once installed and loaded, Seurat can use it automatically.

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
# The marker table is already saved above as
# "data/gbm_top10_markers_per_cluster.csv".
#
# A previous version also wrote a hardcoded absolute Windows path pointing to
# one contributor's Downloads folder, which made the script fail on any other
# machine. That redundant copy has been removed so the script stays portable.
#-------------------------------------------------------------------------------


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
  "LST1", "TYROBP", "AIF1", "PTPRC", "C1QA",
  
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


#===============================================================================
# 4. Broad cell-type annotation
#===============================================================================


#-------------------------------------------------------------------------------
# Start with all clusters as candidate neoplastic
#-------------------------------------------------------------------------------

gbm$broad_celltype <- "Candidate neoplastic"


#-------------------------------------------------------------------------------
# Myeloid
#-------------------------------------------------------------------------------
#
# Cluster 1:
#   TMEM119, FOLR2, C1QA, C1QB, C1QC and CX3CR1,
#   together with canonical LST1, TYROBP, AIF1 and PTPRC,
#   strongly support a microglial/myeloid identity.
#
# Cluster 12:
#   LYZ, S100A8, S100A9, C15orf48 and ADAM8,
#   together with canonical myeloid markers,
#   support an inflammatory myeloid / monocyte-macrophage identity.
#
# Cluster 14:
#   LST1, TYROBP, AIF1, PTPRC and C1QA support myeloid lineage.
#   Strong MKI67, TOP2A, RRM2, TK1 and CEP55 expression indicates
#   that this population is proliferating/cycling.
#
# Therefore, cluster 14 is retained within the broad Myeloid category.

gbm$broad_celltype[
  gbm$seurat_clusters %in% c("1", "12", "14")
] <- "Myeloid"


#-------------------------------------------------------------------------------
# Oligodendrocytes
#-------------------------------------------------------------------------------
#
# Cluster 10:
#   MAG, MOG, OPALIN, NKX6-2, CNDP1 and TMEM125,
#   together with canonical MAG, MOG, MOBP and CLDN11 expression,
#   strongly support mature oligodendrocyte identity.

gbm$broad_celltype[
  gbm$seurat_clusters == "10"
] <- "Oligodendrocytes"


#-------------------------------------------------------------------------------
# T cells
#-------------------------------------------------------------------------------
#
# Cluster 13:
#   CD3D, CD3E, CD3G, CD8B, TRBC1 and GZMA
#   strongly support T-cell identity.

gbm$broad_celltype[
  gbm$seurat_clusters == "13"
] <- "T cells"


#-------------------------------------------------------------------------------
# Fibroblast / stromal
#-------------------------------------------------------------------------------
#
# Cluster 16:
#   COL1A1, COL3A1, LUM, DCN and COL6A3
#   strongly support a fibroblast/stromal identity.
#
# RGS5 and other perivascular features are also present,
# but "Fibroblast" is retained as the broad label.

gbm$broad_celltype[
  gbm$seurat_clusters == "16"
] <- "Fibroblast"


#-------------------------------------------------------------------------------
# Endothelial
#-------------------------------------------------------------------------------
#
# Cluster 17:
#   VWF, SOX17, ECSCR, PCAT19 and MYCT1,
#   together with canonical PECAM1, VWF and CLDN5 expression,
#   strongly support endothelial identity.

gbm$broad_celltype[
  gbm$seurat_clusters == "17"
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
# Final broad cluster mapping:
#
#   Clusters 2, 11, 14  -> Myeloid
#   Cluster 9           -> Oligodendrocytes
#   Cluster 13          -> T cells
#   Cluster 15          -> Fibroblast
#   Cluster 16          -> Endothelial
#
# Remaining clusters are retained as "Candidate neoplastic"
# until further neoplastic-cell validation and state analysis.
#
# Key saved outputs:
#
#   data/gbm_all_cluster_markers.csv
#   data/gbm_top10_markers_per_cluster.csv
#   data/gbm_canonical_marker_dotplot.png
#   data/gbm_broad_celltype_annotation.png
#   data/gbm_broad_celltype_annotated.rds
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