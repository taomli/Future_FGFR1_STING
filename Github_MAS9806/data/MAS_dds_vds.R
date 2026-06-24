library(purrr)
library(ggplot2)
library(pheatmap)
library(DESeq2)

figure.dir <- "/Volumes/T9/Github_MAS9806/MAS9806_RNAseq_dds/QC_plots"

# Read load count data
data <- read.table("mas9806_counts.tsv",
                   sep = "\t",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

# Make Ensembl_ID rownames
rownames(data) = data[,1]

# Remove Ensembl_ID column
count_data = data[,2:31]

# Import colData for all samples (meta)

meta <- read.csv("coldata_MAS98.06.csv", sep=";")
rownames(meta) <- meta$X
meta <- meta[,!grepl("X",colnames(meta))]

# Check that sample names match in both files

all(colnames(count_data) %in% rownames(meta) & rownames(meta) %in% colnames(count_data))
all(colnames(count_data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta, design = ~ treatment)

## Keep genes with ≥10 counts in at least 5 samples

keep <- rowSums(counts(dds) >= 10) >= 5 

dds <- dds[keep,]

dds <- DESeq(dds)

saveRDS(dds, "MAS98.06_dds.RDS")

### Perform variance stabilizing transformation

vsd <- vst(dds, blind=FALSE)

saveRDS(vsd, "MAS98.06_vsd.RDS")


######################################################################
## 1.3 Generate Quality Control Plots and Save to Directory         ##
######################################################################


# Perform rlog transformation for PCA of raw count data
rld <- rlog(dds, blind = FALSE)

# PCA plot for rlog transformed data (raw counts)
pcaData_rlog <- plotPCA(rld, intgroup = "treatment", returnData = TRUE)
percentVar_rlog <- round(100 * attr(pcaData_rlog, "percentVar"))

p1 <- ggplot(pcaData_rlog, aes(x = PC1, y = PC2, color = treatment, label = rownames(colData(dds)))) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.5) +
  xlab(paste0("PC1: ", percentVar_rlog[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_rlog[2], "% variance")) +
  ggtitle("PCA of rlog Transformed Data (Raw Counts)")

ggsave(file.path("PCA_rlog_counts.png"), plot = p1, width = 8, height = 6)

# PCA plot for vst transformed data
pcaData_vst <- plotPCA(vsd, intgroup = "treatment", returnData = TRUE)
percentVar_vst <- round(100 * attr(pcaData_vst, "percentVar"))

p2 <- ggplot(pcaData_vst, aes(x = PC1, y = PC2, color = treatment, label = rownames(colData(dds)))) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.5) +
  xlab(paste0("PC1: ", percentVar_vst[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_vst[2], "% variance")) +
  ggtitle("PCA of VST Transformed Data")

ggsave(file.path(figure.dir, "PCA_vst_counts.png"), plot = p2, width = 8, height = 6)

# Heatmap of sample-to-sample distances (rlog counts)
sampleDists_rlog <- dist(t(assay(rld)))
sampleDistMatrix_rlog <- as.matrix(sampleDists_rlog)
rownames(sampleDistMatrix_rlog) <- paste(rld$treatment, colnames(rld), sep = "-")
colnames(sampleDistMatrix_rlog) <- rownames(sampleDistMatrix_rlog)
heatmap_rlog <- pheatmap(sampleDistMatrix_rlog,
                         clustering_distance_rows = sampleDists_rlog,
                         clustering_distance_cols = sampleDists_rlog,
                         main = "Sample-to-sample distances (rlog counts)")

ggsave(file.path(figure.dir, "heatmap_rlog_counts.png"), plot = heatmap_rlog, width = 8, height = 6)

# Heatmap of sample-to-sample distances (vst counts)
sampleDists_vst <- dist(t(assay(vsd)))
sampleDistMatrix_vst <- as.matrix(sampleDists_vst)
rownames(sampleDistMatrix_vst) <- paste(vsd$treatment, colnames(vsd), sep = "-")
colnames(sampleDistMatrix_vst) <- rownames(sampleDistMatrix_vst)
heatmap_vst <- pheatmap(sampleDistMatrix_vst,
                        clustering_distance_rows = sampleDists_vst,
                        clustering_distance_cols = sampleDists_vst,
                        main = "Sample-to-sample distances (vst counts)")

ggsave(file.path(figure.dir, "heatmap_vst_counts.png"), plot = heatmap_vst, width = 8, height = 6)

# Plot mean-variance relationship before and after VST
meanSdPlot <- function(counts) {
  mns <- rowMeans(counts)
  vars <- rowVars(counts)
  data.frame(mean = mns, variance = vars)
}

raw_plot_data <- meanSdPlot(counts(dds))
vst_plot_data <- meanSdPlot(assay(vsd))

p3 <- ggplot(raw_plot_data, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.3) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Mean-Variance Relationship (Raw Counts)")

ggsave(file.path(figure.dir, "mean_variance_raw_counts.png"), plot = p3, width = 8, height = 6)

p4 <- ggplot(vst_plot_data, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.3) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Mean-Variance Relationship (VST Counts)")

ggsave(file.path(figure.dir, "mean_variance_vst_counts.png"), plot = p4, width = 8, height = 6)

