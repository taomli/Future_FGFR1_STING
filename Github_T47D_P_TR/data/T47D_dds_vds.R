library(DESeq2)
library(ggplot2)
library(pheatmap)

########### Set seed ##############

set.seed(42)

########################################################
######################## 1.0 Read data #################
########################################################

# Define various directories 

data.dir = ("/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/")

setwd(data.dir)

count_data <- read.table(
  "/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/T47D_counts.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Import colData for all samples (meta)

coldata <- read.csv("/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/T47D_coldata.csv", sep=";", row.names = 1)

# Reorder coldata according to colnames of count_data
coldata <- coldata[colnames(count_data), , drop = FALSE]
all(rownames(coldata) == colnames(count_data))

#### Raw count data matrix for TR01

count_data_treatment <- count_data[,!grepl("PV",colnames(count_data))]

#### Create coldata for TR01 samples (meta_tr)

coldata_treatment <- coldata[!grepl("PV",rownames(coldata)),]

#### Raw count data matrix for PAR/TR comparison

count_data_cell_line <- count_data[,c("PV1","PV2","PV3","TV1","TV2","TV3")]

#### Create coldata for TR01 samples (meta_tr)

coldata_cell_line <- coldata[c("PV1","PV2","PV3","TV1","TV2","TV3"),]



#################################################################################################
## Create DESeq2Dataset (dds, and vsd) object for parental and resistant vehicle samples ####
#################################################################################################

dds_cell_line <- DESeqDataSetFromMatrix(countData = count_data_cell_line, colData = coldata_cell_line, design = ~ cell_line)

keep <- rowSums(counts(dds_cell_line) >= 10) >= 3 ## Remove rows with less than 10 reads total

dds_cell_line <- dds_cell_line[keep,]

dds_cell_line <- DESeq(dds_cell_line)

saveRDS(dds_cell_line, "T47D_dds_cell_line.RDS")

### Perform variance stabilizing transformation

vsd_cell_line <- vst(dds_cell_line, blind=FALSE)

saveRDS(vsd_cell_line, "T47D_vsd_cell_line.RDS")


######################################################################################################
######### 1.3 Create DESeq2Dataset (dds, and vsd) object for treatments of resistant cells ###########
######################################################################################################

dds_treatment <- DESeqDataSetFromMatrix(countData = count_data_treatment, colData = coldata_treatment, design = ~ treatment)

keep <- rowSums(counts(dds_treatment) >= 10) >= 3

dds_treatment <- dds_treatment[keep,]

dds_treatment <- DESeq(dds_treatment)

saveRDS(dds_treatment, "T47D_dds_treatment.RDS")

### Transform counts for data visualization

vsd_treatment <- vst(dds_treatment, blind=FALSE)

saveRDS(vsd_treatment, "T47D_vsd_treatment.RDS")

