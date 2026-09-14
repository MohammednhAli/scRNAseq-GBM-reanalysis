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

# Before running this script, obtain the raw count matrix for
# GEO accession GSE229779 and place it in the local data/ folder:
#
# data/GSE229779_countsMatrix.tsv.gz
#
# The data/ folder is excluded from this GitHub repository via
# .gitignore, so the raw dataset must be downloaded separately
# by each user.
#
# fread() reads the compressed tab-separated count matrix.
#
# Expected structure:
#
# gene      Cell1   Cell2   Cell3 ...
# TP53        0       2       0
# EGFR        5       1       8
# SOX2        2       0       3
#
# Rows    = genes
# Columns = cells
# Values  = raw expression counts

counts_data <- fread("data/GSE229779_countsMatrix.tsv.gz")


# Check the size of the imported table
dim(counts_data)

# View the first 5 rows and first 5 columns
counts_data[1:5, 1:5]


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

# Check the first 5 column names
names(counts_data)[1:5]


# Save the gene names separately before removing the text column.
#
# Example:
# gene_names =
# MIR1302-2HG
# AL627309.1
# AL669831.5

gene_names <- counts_data$gene

# Check the first few gene names
head(gene_names)


# Remove the gene-name column from the count table.
#
# We do this because the expression matrix itself should contain
# only numerical count values.

counts_data[, gene := NULL]

# Check the first 5 remaining columns
head(counts_data[, 1:5])


# Convert the data.table into a standard R matrix.
#
# Seurat requires a gene x cell count matrix as input.

counts_matrix <- as.matrix(counts_data)


# Add the gene names back as row names.
rownames(counts_matrix) <- gene_names

# The matrix now has:
#   rows    = genes
#   columns = cells
#   values  = raw expression counts
#
# Example structure:
#
#             GBM21_AAACCT...   GBM21_AAACCT...   GBM21_AAACCT...
# MIR1302-2HG        0                 0                 0
# AL627309.1         0                 0                 0
# AL669831.5         0                 0                 0
# FAM87B             0                 0                 0
# LINC00115          1                 0                 0

# View the first 5 genes and first 5 cells
counts_matrix[1:5, 1:5]



# ============================================================
# 5. CHECK THE COUNT MATRIX
# ============================================================

# Check number of genes x number of cells
dim(counts_matrix)

# Check the first 5 gene names
rownames(counts_matrix)[1:5]

# Check the first 5 cell names
colnames(counts_matrix)[1:5]

# Confirm that the object is a matrix
class(counts_matrix)


# ============================================================
# 6. REMOVE TEMPORARY OBJECT
# ============================================================

# Remove counts_data because counts_matrix is now the object
# used for downstream analysis
rm(counts_data)

# Free unused memory
gc()

# ============================================================
# 7. CREATE RAW SEURAT OBJECT
# ============================================================

# Create a Seurat object from the raw gene x cell count matrix.
#
# min.cells = 0 and min.features = 0 prevent filtering during
# object creation, so all genes and cells are retained.
#
# QC and filtering will be performed in the next analysis script.

gbm <- CreateSeuratObject(
  counts = counts_matrix,
  project = "GBM",
  min.cells = 0,
  min.features = 0
)

# ============================================================
# 8. CHECK THE SEURAT OBJECT
# ============================================================

# Print a summary of the Seurat object
gbm

# Check number of genes x number of cells
dim(gbm)

# ============================================================
# 9. CLEAN UP MEMORY
# ============================================================

# Remove temporary objects that are no longer needed
rm(counts_matrix, gene_names)

# Free unused memory
gc()

# ============================================================
# 10. SAVE RAW SEURAT OBJECT
# ============================================================

# Save the unfiltered Seurat object as an RDS file.
#
# This avoids re-reading the large count matrix and rebuilding
# the Seurat object every time the next analysis script is run.
#
# The saved object contains the raw counts and will be used as
# the input for QC and doublet detection.

saveRDS(
  gbm,
  file = "data/gbm_raw_seurat.rds"
)
