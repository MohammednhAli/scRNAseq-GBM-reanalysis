# ============================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
#
# Script: 01_QC_DoubletDetection.R
#
# Purpose:
#   1. Calculate and visualize QC metrics
#   2. Filter low-quality cells
#   3. Detect doublets using scDblFinder
#   4. Inspect and remove doublets
# ============================================================


# ============================================================
# 1. INSTALL PACKAGES IF NEEDED
# ============================================================

# Uncomment these lines only if the packages are not installed

# install.packages("ggplot2")
# install.packages("BiocManager")

# BiocManager::install("SingleCellExperiment")
# BiocManager::install("scDblFinder")


# ============================================================
# 2. LOAD RAW SEURAT OBJECT
# ============================================================

# Load the raw Seurat object created by 00_LoadData.R.
#
# readRDS() loads the saved Seurat object directly, avoiding
# the need to re-read the large count matrix and recreate the
# Seurat object every time the QC script is run.

gbm <- readRDS(
  "data/gbm_raw_seurat.rds"
)


# ============================================================
# 3. LOAD PACKAGES
# ============================================================

library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)


# ============================================================
# 4. QUALITY CONTROL
# ============================================================

# Calculate the percentage of mitochondrial RNA in each cell.
#
# PercentageFeatureSet() finds genes matching the pattern
# "^MT-" and calculates what percentage of the total RNA
# counts comes from mitochondrial genes.

gbm[["percent.mt"]] <- PercentageFeatureSet(
  gbm,
  pattern = "^MT-"
)


# ============================================================
# 5. VISUALIZE QC METRICS
# ============================================================

# VlnPlot() shows the distribution of the three main
# QC measurements:
#
# nFeature_RNA = number of detected genes per cell
#
# nCount_RNA = total RNA counts per cell
#
# percent.mt = percentage of mitochondrial RNA
#
# These measurements are examined together to identify
# low-quality or unusual cells.

VlnPlot(
  gbm,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)


# ------------------------------------------------------------
# Total RNA counts vs detected genes
# ------------------------------------------------------------

# FeatureScatter() plots two QC measurements against each
# other.
#
# This allows us to see the relationship between the total
# amount of RNA in a cell and the number of genes detected.

FeatureScatter(
  gbm,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)


# ------------------------------------------------------------
# Total RNA counts vs mitochondrial RNA
# ------------------------------------------------------------

# This plot shows the relationship between total RNA counts
# and mitochondrial RNA percentage.

FeatureScatter(
  gbm,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)


# ============================================================
# 6. CHECK QC RANGES
# ============================================================

# Check the observed range of detected genes per cell.

range(gbm$nFeature_RNA)


# Check the observed range of mitochondrial RNA percentage.

range(gbm$percent.mt)


# Observed:
#
# nFeature_RNA = 212 - 4499
# percent.mt   = 0 - 12.49851
#
# The QC violin plots show clear upper and lower boundaries,
# with no cells extending beyond these ranges.
#
# Together with the observed minimum and maximum values, this
# is consistent with the GEO count matrix having already
# undergone QC filtering before deposition.
#
# Therefore, the filtering step below mainly confirms the
# existing QC boundaries rather than removing additional cells.


# ============================================================
# 7. FILTER LOW-QUALITY CELLS
# ============================================================

# subset() keeps a cell only when every condition is TRUE.
#
# Here we keep cells with:
#
# more than 200 detected genes
# fewer than 4500 detected genes
# less than 12.5% mitochondrial RNA
#
# The "&" symbol means AND, so all three conditions must
# be satisfied for a cell to remain.
#
# The deposited GSE229779 count matrix had already undergone
# QC filtering in the original study.
#
# In our dataset, the observed ranges were:
#
# nFeature_RNA = 212 - 4499
# percent.mt   = 0 - 12.49851
#
# Therefore, all cells already satisfy the criteria below.
# The filtering command is retained here to demonstrate the
# QC step and confirm the existing filtering boundaries.

gbm <- subset(
  gbm,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 4500 &
    percent.mt < 12.5
)


# Check the dimensions after QC filtering.
#
# Rows    = genes
# Columns = cells

dim(gbm)


# ============================================================
# 8. CONVERT TO SingleCellExperiment
# ============================================================

# scDblFinder works with a SingleCellExperiment object.
#
# as.SingleCellExperiment() converts the Seurat object into
# the format required by scDblFinder.

sce_gbm <- as.SingleCellExperiment(gbm)


# ============================================================
# 9. DETECT DOUBLETS
# ============================================================

# A doublet occurs when two cells are captured in the same
# droplet and are measured as if they were one cell.
#
# scDblFinder() scores each cell and predicts whether it is
# a singlet or a doublet.


# set.seed() makes the result reproducible.

set.seed(100)


# If scDblFinder is run without specifying the sample identity,
# it will treat all cells as if they came from one combined sample.
#
# This is not appropriate here because cells from different samples
# were processed separately and could not have formed doublets
# with cells from another sample.
#
# Therefore, doublet detection is performed separately within
# each sample to avoid impossible cross-sample doublets and
# reduce incorrect classifications.

sce_gbm <- scDblFinder(
  sce_gbm,
  samples = "orig.ident"
)


# ============================================================
# 10. ADD DOUBLET RESULTS BACK TO SEURAT
# ============================================================

# scDblFinder produces two important outputs:
#
# scDblFinder.score
#   = a continuous doublet score
#
# scDblFinder.class
#   = classification as "singlet" or "doublet"
#
# Add both results to the Seurat metadata.

gbm$doublet_score <-
  colData(sce_gbm)$scDblFinder.score

gbm$doublet_class <-
  colData(sce_gbm)$scDblFinder.class


# ============================================================
# 11. INSPECT DOUBLET RESULTS
# ============================================================

# Count the number of predicted singlets and doublets.

table(gbm$doublet_class)


# Check the doublet classification within each sample.

table(
  gbm$orig.ident,
  gbm$doublet_class
)


# ============================================================
# 12. VISUALIZE DOUBLET DETECTION
# ============================================================

# ------------------------------------------------------------
# RNA counts vs doublet score
# ------------------------------------------------------------

# Plot nCount_RNA against the doublet score.
#
# Cells predicted as doublets are shown separately from
# singlets.
#
# facet_wrap() shows the result separately for each sample.

ggplot(
  gbm@meta.data,
  aes(
    x = nCount_RNA,
    y = doublet_score,
    color = doublet_class
  )
) +
  geom_point(
    size = 0.7,
    alpha = 0.5
  ) +
  facet_wrap(~ orig.ident) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score"
  ) +
  theme_classic()


# ------------------------------------------------------------
# Compare doublet scores
# ------------------------------------------------------------

# VlnPlot() compares the doublet-score distributions between
# cells classified as singlets and cells classified as
# doublets.

VlnPlot(
  gbm,
  features = "doublet_score",
  group.by = "doublet_class",
  split.by = "orig.ident"
)


# ============================================================
# 13. REMOVE DOUBLETS
# ============================================================

# Keep only cells classified as singlets.

gbm <- subset(
  gbm,
  subset = doublet_class == "singlet"
)


# Check the dimensions after removing doublets.

dim(gbm)


# Confirm that only singlets remain.

table(gbm$doublet_class)


# Remove the unused "doublet" factor level.

gbm$doublet_class <- droplevels(
  gbm$doublet_class
)


# ============================================================
# 14. VISUALIZE AFTER DOUBLET REMOVAL
# ============================================================

# All remaining cells are singlets.
#
# Compare the remaining doublet-score distributions across
# samples.

VlnPlot(
  gbm,
  features = "doublet_score",
  group.by = "orig.ident"
)


# ============================================================
# 15. FINAL QC-PASSED OBJECT
# ============================================================

# gbm now contains:
#
# QC-passed cells
#        +
# predicted singlets
#
# This object can now be used for the next stage of the
# scRNA-seq analysis.

gbm


# Check the final number of genes and cells.

dim(gbm)


# ============================================================
# 16. SAVE FINAL QC-PASSED OBJECT
# ============================================================

# Save the final QC-passed Seurat object as an RDS file.
#
# saveRDS() stores the complete R object on disk, including:
#
#   - raw RNA counts
#   - QC metadata
#   - doublet scores
#   - retained singlet cells
#
# This allows downstream scripts to load the completed QC object
# directly without repeating data loading, QC calculations or
# scDblFinder doublet detection.
#
# The saved object will be used as the input for normalization.

saveRDS(
  gbm,
  file = "data/gbm_QC_singlets.rds"
)
