########
#
# Analysis of Bulk RNA-Seq Data using DESeq2
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
#
########

################################################################################
#### PART 1: Introduction to DESeq2 ####
################################################################################

################################################################################
#### SECTION 1: Loading Packages ####
################################################################################

# DESeq2 is one of the most popular packages for differential expression analysis.
# It handles normalization, dispersion estimation, and statistical testing
# all in one streamlined workflow.

library(DESeq2)
library(tidyverse)
library(tidymodels)

################################################################################
#### SECTION 2: Loading and Preparing the Data ####
################################################################################

# Load the same TCGA expression data we used in 3.1.1:
d <- readRDS("Data/TCGA_Formatted_data.rds")

# We select two cancer types to compare (Breast Cancer vs. Colorectal Adenocarcinoma):
d <- filter(d, Cancertype %in% c("BRCA", "COAD"))

# Clean the Sample ID (remove extra characters):
d <- mutate(d, SampleID = str_sub(Sample, end = -14))

# Store the cancer type information before removing it:
sample_info <- d %>%
dplyr::select(SampleID, Cancertype)

# Keep only SampleID and gene expression columns:
d <- dplyr::select(d, SampleID, everything()) %>%
dplyr::select(-c(Sample, Cancertype))

################################################################################
#### SECTION 3: Filtering Low-Expression Genes ####
################################################################################

# Reformat to long format and count the number of samples with expression > 0:
Summary_d <- d %>%
pivot_longer(cols = -SampleID) %>%
group_by(name) %>%
mutate(n_above_0 = sum(value > 0)) %>%
distinct(name, n_above_0) %>%
dplyr::rename(gene_id = name) %>%
ungroup()

# We keep only genes that are expressed in at least 25% of samples:
genes_to_include <- filter(Summary_d, n_above_0 > 0.25 * ncol(d))$gene_id

# Subset the data to include only these genes:
d <- dplyr::select(d, SampleID, any_of(genes_to_include))

# EXERCISE A: Remaining Genes
# How many genes remain after filtering?
# Write your answer here:
#
#

################################################################################
#### SECTION 4: Selecting Most Variable Genes ####
################################################################################

# To reduce computation time, we keep only the 200 most variable genes.
# We calculate variance across all samples for each gene, then select the top 200.

# Calculate variance for each gene (column) in the expression matrix:
gene_variances <- d %>%
dplyr::select(-SampleID) %>%
summarise(across(everything(), \(x) var(x, na.rm = TRUE))) %>%
pivot_longer(cols = everything(), names_to = "gene", values_to = "variance")

# Select top 200 most variable genes:
top_variable_genes <- gene_variances %>%
slice_max(order_by = variance, n = 200) %>%
pull(gene)

# Subset data to keep only the most variable genes:
d <- dplyr::select(d, SampleID, all_of(top_variable_genes))

# STATISTICAL NOTE:
# This approach is statistically valid ("unsupervised" feature selection) because
# we select genes based on overall variance across ALL samples, without using
# the cancer type labels (BRCA vs COAD). This means:
#   1. We are NOT peeking at the outcome variable when selecting features.
#   2. The selection is based purely on gene expression variability.
#   3. This avoids "data leakage" that would occur if we selected genes based
#      on their differential expression between groups.

################################################################################
#### SECTION 5: Preparing DESeq2 Input ####
################################################################################

# DESeq2 requires three pieces of information:
# 1. countData: Raw counts in matrix format (genes as rows, samples as columns)
# 2. colData: Metadata about the samples (columns) in countData
# 3. design: A formula specifying which variable(s) to test

# Remove duplicate SampleIDs (keep only the first occurrence):
n_before <- nrow(d)
d <- d %>%
  distinct(SampleID, .keep_all = TRUE)
n_after <- nrow(d)
cat("Removed", n_before - n_after, "duplicate samples. Remaining:", n_after, "\n")

sample_info <- sample_info %>%
  distinct(SampleID, .keep_all = TRUE)

# Format the count matrix (DESeq2 expects genes as rows, samples as columns):
count_matrix <- d %>%
column_to_rownames("SampleID") %>%
t() %>%  # Transpose: genes become rows, samples become columns
as.matrix()

# Check dimensions:
dim(count_matrix)  # Should be: genes x samples

# Format the colData (sample metadata):
colData <- sample_info %>%
filter(SampleID %in% colnames(count_matrix)) %>%
column_to_rownames("SampleID")

# Make sure the order of rows in colData matches the order of columns in count_matrix:
colData <- colData[colnames(count_matrix), , drop = FALSE]

# Convert Cancertype to factor:
colData$Cancertype <- as.factor(colData$Cancertype)

# EXERCISE B: Understanding the Data Structure
# Why does DESeq2 require genes as rows and samples as columns?
# Why is it important that colData rows match count_matrix columns?
# Write your answer here:
#
#

################################################################################
#### SECTION 6: Creating the DESeq2 Dataset ####
################################################################################

# Create the DESeqDataSet object:
dds <- DESeqDataSetFromMatrix(
countData = count_matrix,
colData   = colData,
design    = ~ Cancertype
)

# Inspect the object:
dds

# CHECKPOINT: The design formula "~ Cancertype" tells DESeq2 to test for
# differences between cancer types (BRCA vs COAD).

################################################################################
#### SECTION 7: Running DESeq2 ####
################################################################################

# Run the DESeq2 pipeline (this performs normalization, dispersion estimation,
# and statistical testing all in one step):
dds <- DESeq(dds)

# Extract results:
res <- results(dds)
res

# EXERCISE C: Understanding the Output
# What do the columns in the results mean?
# - baseMean: ?
# - log2FoldChange: ?
# - lfcSE: ?
# - stat: ?
# - pvalue: ?
# - padj: ?
# Write your answers here:
#
#

################################################################################
#### SECTION 8: Exploring the Results ####
################################################################################

# Convert to data frame for easier manipulation:
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df %>%
dplyr::select(gene_id, everything()) %>%
arrange(padj)

head(res_df)

# How many genes are significantly differentially expressed (padj < 0.05)?
sum(res_df$padj < 0.05, na.rm = TRUE)

# EXERCISE D: Top Differentially Expressed Gene
# Which gene is the most significantly differentially expressed?
# Is it up- or down-regulated in BRCA compared to COAD?
# Write your answer here:
#
#

################################################################################
#### SECTION 9: Volcano Plot ####
################################################################################

# Create a volcano plot to visualize the results:
# Replace padj = 0 with a small value to avoid -log10(0) = Inf
res_df <- res_df %>%
  mutate(
    padj_plot = ifelse(padj == 0 | is.na(padj), .Machine$double.xmin, padj),
    neg_log10_padj = -log10(padj_plot)
  )

# Identify the top 30 most significant genes for labeling:
top_genes <- res_df %>%
  slice_min(order_by = padj, n = 30)

ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
geom_point(aes(color = padj < 0.05), alpha = 0.5) +
scale_color_manual(values = c("grey50", "red"), labels = c("Not Sig.", "FDR < 0.05")) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
geom_text(data = top_genes, aes(label = gene_id), 
          size = 3, hjust = -0.1, vjust = 0.5, check_overlap = TRUE) +
labs(
  title = "Volcano Plot: DESeq2 Differential Expression (BRCA vs COAD)",
  x = "Log2 Fold-Change",
  y = "-log10(adjusted p-value)",
  color = "Significance"
) +
theme_minimal(base_size = 30)

# EXERCISE E: Comparing to Manual Analysis
# Compare this volcano plot to the one from 3.1.1.
# Do you see similar patterns? Are the top genes the same?
# Write your answer here:
#
#

################################################################################
#### SECTION 10: MA Plot ####
################################################################################

# DESeq2 also provides an MA plot, which shows log fold-change vs. mean expression:
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)  # Increase font sizes
plotMA(res, ylim = c(-5, 5))

# EXERCISE F: Interpreting the MA Plot
# What does this plot show? Why might genes with low expression have
# more variable fold-changes?
# Write your answer here:
#
#

################################################################################
#### KEY TAKEAWAYS ####
################################################################################

# 1. DESeq2 provides a complete pipeline for differential expression analysis,
#    handling normalization and statistical testing automatically.
#
# 2. DESeq2 uses a negative binomial model, which is appropriate for count data
#    with overdispersion (variance > mean).
#
# 3. The adjusted p-value (padj) accounts for multiple testing using the
#    Benjamini-Hochberg procedure (FDR).
#
# 4. Volcano plots and MA plots are useful for visualizing differential
#    expression results.
#
# 5. Comparing results from manual analysis (3.1.1) and DESeq2 helps validate
#    findings and understand the strengths of each approach.



