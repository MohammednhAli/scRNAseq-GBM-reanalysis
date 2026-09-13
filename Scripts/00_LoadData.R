# ============================================================
# GBM scRNA-seq ANALYSIS
# Dataset: GSE229779
# Drug and single-cell gene expression integration identifies
# sensitive and resistant glioblastoma cell populations
#
# Script: 00_LoadData.R
# Purpose:
#   1. Load the raw count matrix
#   2. Convert it into a gene x cell matrix
#   3. Create an unfiltered Seurat object
#
# No QC or filtering is performed in this script.
# ============================================================


# ============================================================
# 1. INSTALL PACKAGES IF NEEDED
# ============================================================

# Uncomment these lines only if the packages are not installed
# install.packages("data.table")
# install.packages("Seurat")


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

# data.table provides:
#   fread()    -> fast reading of large tabular files
#   setnames() -> simple renaming of columns
#
# It is useful here because the scRNA-seq count matrix is large
# and fread() is faster and more memory-efficient than read.table().

library(data.table)

# Seurat is used to create the scRNA-seq object that will later
# be used for QC, normalization, clustering and annotation.

library(Seurat)


# ============================================================
# 3. READ THE COUNT MATRIX
# ============================================================

# fread() reads the compressed tab-separated count matrix.
#
# Expected structure:
#
# gene      Cell1   Cell2   Cell3 ...
# TP53        0       2       0
# EGFR        5       1       8
# SOX2        2       0       3
#
# Rows = genes
# Columns = cells
# Values = raw expression counts

counts_data <- fread("data/GSE229779_countsMatrix.tsv.gz")


# Check the size of the imported table
dim(counts_data)

# Look at the first few rows
head(counts_data)


# ============================================================
# 4. PREPARE GENE x CELL COUNT MATRIX
# ============================================================

# The first column contains gene names.
#
# When fread() reads the file, this column may have a generic
# name such as V1.
#
# setnames() renames column 1 to "gene" so that the code is
# easier to understand.

setnames(counts_data, 1, "gene")


# Save the gene names separately before removing the text column.
#
# Example:
# gene_names =
# TP53
# EGFR
# SOX2

gene_names <- counts_data$gene

# Check the first gene names
head(gene_names)


# Remove the gene-name column from the count table.
#
# We do this because the expression matrix itself should contain
# only numerical count values.

counts_data[, gene := NULL]


# Convert the data.table into a standard R matrix.
#
# Seurat requires a gene x cell count matrix as input.

counts_matrix <- as.matrix(counts_data)


# Add the gene names back as row names.
#
# The matrix will now look like:
#
#        Cell1   Cell2   Cell3
# TP53      0       2       0
# EGFR      5       1       8
# SOX2      2       0       3

rownames(counts_matrix) <- gene_names


# ============================================================
# 5. CHECK THE COUNT MATRIX
# ============================================================

# Number of genes x number of cells
dim(counts_matrix)

# Check first gene names
rownames(counts_matrix)[1:5]

# Check first cell names
colnames(counts_matrix)[1:5]

# Confirm that the object is a matrix
class(counts_matrix)


# ============================================================
# 6. REMOVE TEMPORARY OBJECT
# ============================================================

# counts_data is no longer needed because we now have
# counts_matrix.

rm(counts_data)

# Ask R to free unused memory
gc()


# ============================================================
# 7. CREATE RAW SEURAT OBJECT
# ============================================================

# CreateSeuratObject() converts the raw count matrix into a
# Seurat object for downstream single-cell analysis.
#
# min.cells = 0 and min.features = 0 are used because this script
# should NOT perform QC or filtering.
#
# Cell and gene filtering will be handled in the QC script.

gbm <- CreateSeuratObject(
  counts = counts_matrix,
  project = "GBM",
  min.cells = 0,
  min.features = 0
)


# ============================================================
# 8. CHECK THE SEURAT OBJECT
# ============================================================

# Print summary of the Seurat object
gbm

# Check number of genes x cells
dim(gbm)