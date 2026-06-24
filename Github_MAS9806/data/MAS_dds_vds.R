library(DESeq2)

data.dir = "/Volumes/T9/Github_projects/Future_FGFR1_STING/Github_MAS9806/data/"

# Read load count data
data <- read.table(paste0(data.dir,"mas9806_counts.tsv"),
                   sep = "\t",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

rownames(data) <- data$Ensembl_ID
count_data <- data[, !colnames(data) %in% "Ensembl_ID"]

# Import colData for all samples (meta)

meta <- read.csv(file.path(data.dir, "coldata_MAS98.06.csv"), sep=";")
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

saveRDS(dds, file.path(data.dir, "MAS98.06_dds.RDS"))
### Perform variance stabilizing transformation

vsd <- vst(dds, blind=FALSE)

saveRDS(vsd, file.path(data.dir, "MAS98.06_vsd.RDS"))
