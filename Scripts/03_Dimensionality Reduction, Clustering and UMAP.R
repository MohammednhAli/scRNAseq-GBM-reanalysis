#===============================================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
#
# Script: 03_Dimensionality_Reduction_Integration_Clustering_UMAP.R
#
# Purpose:
#   1. Load the normalized/scaled Seurat object.
#   2. Run PCA.
#   3. Select the number of PCs to retain.
#   4. Generate UMAP before integration.
#   5. Perform Harmony integration across samples.
#   6. Construct the neighborhood graph using Harmony.
#   7. Cluster cells.
#   8. Generate UMAP after Harmony integration.
#   9. Compare UMAP before vs after integration.
#  10. Save the final processed Seurat object.
#===============================================================================


#===============================================================================
# 0. Load packages and processed Seurat object
#===============================================================================

library(Seurat)
library(harmony)

gbm <- readRDS("data/gbm_normalized_HVG_scaled.rds")

# Check that the object loaded correctly
gbm


#===============================================================================
# 1. Run PCA
#===============================================================================

# Purpose:
# PCA reduces the dimensionality of the dataset.
#
# Instead of representing each cell using thousands of genes, PCA summarizes
# the major directions of variation using principal components (PCs).
#
# PCA is calculated using the highly variable genes selected in the previous
# normalization/HVG step.

gbm <- RunPCA(
  gbm,
  features = VariableFeatures(object = gbm)
)


#-------------------------------------------------------------------------------
# Visualize PCA
#-------------------------------------------------------------------------------

DimPlot(
  gbm,
  reduction = "pca"
)

# Result:
# The plot shows PC1 versus PC2.
#
# PC1 and PC2 are the first two major directions of variation across cells.
# Each PC is a weighted combination of many genes.
#
# This plot provides an initial view of large-scale structure and
# sample-associated differences in the dataset.


#===============================================================================
# 2. Select the number of PCs
#===============================================================================


#-------------------------------------------------------------------------------
# 2A. Elbow plot
#-------------------------------------------------------------------------------

ElbowPlot(gbm)

# Purpose:
# The elbow plot shows the relative contribution of successive PCs.
#
# Earlier PCs explain more variation, while later PCs contribute progressively
# less.
#
# Based on the elbow plot and inspection of the PC heatmaps, PCs 1-20 were
# retained for downstream analysis.
#
# Selected dimensions:
# dims = 1:20
#
# Important:
# PCs after PC20 are not necessarily noise; they simply explain less variation.


#-------------------------------------------------------------------------------
# 2B. Visualize genes contributing to selected PCs
#-------------------------------------------------------------------------------

DimHeatmap(
  gbm,
  dims = 1:20,
  cells = 500,
  balanced = TRUE
)

# Purpose:
# DimHeatmap visualizes genes with strong positive and negative loadings
# for each selected PC.
#
# dims = 1:20:
#   visualize PCs 1 through 20.
#
# cells = 500:
#   display 500 cells for easier visualization.
#
# balanced = TRUE:
#   display a balanced number of genes with positive and negative loadings.
#
# Result:
# PCs 1-20 showed structured gene-expression patterns across cells, supporting
# their use for downstream integration, neighborhood construction, clustering,
# and UMAP.


#===============================================================================
# 3. UMAP BEFORE integration
#===============================================================================

# Purpose:
# Generate a diagnostic UMAP using the original PCA representation before
# Harmony correction.
#
# This allows direct visualization of sample-associated structure before
# integration.

gbm <- RunUMAP(
  gbm,
  reduction = "pca",
  dims = 1:20,
  reduction.name = "umap_before",
  reduction.key = "UMAPbefore_"
)


#-------------------------------------------------------------------------------
# Visualize samples before integration
#-------------------------------------------------------------------------------

p_before <- DimPlot(
  gbm,
  reduction = "umap_before",
  group.by = "orig.ident"
) +
  ggtitle("Before Harmony integration")

p_before

# Result:
# Before integration, several cell populations were strongly associated with
# individual GBM samples.
#
# This indicates substantial sample-associated structure.
#
# However, sample-specific separation does not automatically indicate technical
# batch effect because genuine biological differences between GBM patients can
# also contribute.


#===============================================================================
# 4. Harmony integration
#===============================================================================

# Purpose:
# Harmony corrects the PCA representation for sample-associated variation.
#
# orig.ident identifies the six individual GBM samples:
#
# GBM21
# GBM41
# GBM47
# GBM49
# GBM51
# GBM53
#
# Harmony uses the PCA coordinates and creates a new corrected dimensional
# representation called "harmony".

gbm <- RunHarmony(
  gbm,
  "orig.ident"
)


# Result:
# A new reduction called "harmony" is stored inside the Seurat object.
#
# The original PCA is still preserved.
#
# Main reductions now include:
#
# pca
# umap_before
# harmony


#===============================================================================
# 5. Construct the neighborhood graph
#===============================================================================

# Purpose:
# Build a nearest-neighbor graph using the Harmony-corrected representation.
#
# Cells that are transcriptionally similar are connected in this graph.
#
# The same selected dimensions, 1-20, are used.

gbm <- FindNeighbors(
  gbm,
  reduction = "harmony",
  dims = 1:20
)

# Result:
# Seurat stores the KNN/SNN graphs in the object.
#
# These graphs are used in the next step for graph-based clustering.


#===============================================================================
# 6. Cell clustering
#===============================================================================

# Purpose:
# Group transcriptionally similar cells using the neighborhood graph.
#
# resolution controls clustering granularity:
#
# lower resolution  -> fewer, broader clusters
# higher resolution -> more, smaller clusters
#
# A moderate resolution of 0.5 is used here.

gbm <- FindClusters(
  gbm,
  resolution = 0.5
)

# Result:
# Each cell receives a cluster identity.
#
# Cluster assignments are stored in:
#
# gbm$seurat_clusters


#===============================================================================
# 7. UMAP AFTER Harmony integration
#===============================================================================

# Purpose:
# Generate the final UMAP using the Harmony-corrected dimensional
# representation.
#
# This UMAP is used for visualization of integrated samples and clusters.

gbm <- RunUMAP(
  gbm,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap_after",
  reduction.key = "UMAPafter_"
)


#-------------------------------------------------------------------------------
# Visualize samples after Harmony integration
#-------------------------------------------------------------------------------

p_after <- DimPlot(
  gbm,
  reduction = "umap_after",
  group.by = "orig.ident"
) +
  ggtitle("After Harmony integration")

p_after

# Result:
# Harmony improves mixing between several samples while preserving distinct
# biological populations.


#-------------------------------------------------------------------------------
# Visualize Seurat clusters
#-------------------------------------------------------------------------------

DimPlot(
  gbm,
  reduction = "umap_after",
  label = TRUE
)

# Result:
# Cells are displayed in two-dimensional UMAP space and labeled according to
# their Seurat cluster identity.


#===============================================================================
# 8. Compare BEFORE vs AFTER Harmony integration
#===============================================================================

p_before + p_after

# Purpose:
# Directly compare sample-associated structure before and after Harmony.
#
# BEFORE:
#
# PCA
#  |
#  v
# UMAP
#
#
# AFTER:
#
# PCA
#  |
#  v
# Harmony
#  |
#  v
# UMAP
#
#
# Interpretation:
# Harmony reduces substantial sample-associated structure.
#
# Complete mixing of all six samples is not required because GBM contains real
# inter-patient biological heterogeneity.


#===============================================================================
# 9. Optional checks
#===============================================================================

# Check stored dimensional reductions
Reductions(gbm)

# Check number of cells in each cluster
table(gbm$seurat_clusters)

# Check sample contribution to each cluster
table(
  gbm$seurat_clusters,
  gbm$orig.ident
)


#===============================================================================
# 10. Save final processed Seurat object
#===============================================================================

saveRDS(
  gbm,
  file = "data/gbm_harmony_clustered.rds"
)

# Result:
# The saved object contains:
#
# - normalized RNA data
# - scaled RNA data
# - highly variable genes
# - PCA
# - UMAP before integration
# - Harmony corrected reduction
# - neighborhood graph
# - Seurat cluster assignments
# - UMAP after Harmony integration

