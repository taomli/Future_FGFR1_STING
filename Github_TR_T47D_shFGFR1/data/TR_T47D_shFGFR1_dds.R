library(DESeq2)

# Read data
counts <- read.delim("TR_T47D_shFGFR1_gene_count.txt", sep = "\t", header = TRUE)

# Extract count matrix
count_matrix <- counts[, paste0("TR_seq_", 1:8)]
rownames(count_matrix) <- counts$gene_id

# Read and format colData
colData <- read.csv2("TR_T47D_shFGFR1_coldata.csv", row.names = 1)
colData$treatment <- factor(colData$treatment, levels = c("sh_CTR_h2o", "sh_FGFR1_h2o"))
colData <- colData[, "treatment", drop = FALSE]

# Verify sample order
stopifnot(all(rownames(colData) == colnames(count_matrix)))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ treatment)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Save objects
saveRDS(dds, file = "TR_T47D_shFGFR1_dds.rds")
