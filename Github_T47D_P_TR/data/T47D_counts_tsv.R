library(purrr)
library(dplyr)

data.dir = ("/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/30_featureCounts/")

# Check that ensemble IDs are equal across all 30 files

setwd(data.dir)
files <- list.files(pattern = "hisat2_FC$")

# Verify that geneids match in all FC files
ref_geneid <- read.table(
  files[1],
  sep = "\t",
  skip = 1,
  header = TRUE,
  stringsAsFactors = FALSE
)$Geneid

all(sapply(files, function(f) {
  all(
    read.table(
      f,
      sep = "\t",
      skip = 1,
      header = TRUE,
      stringsAsFactors = FALSE
    )$Geneid == ref_geneid
  )
}))

# List colnames of FC file
colnames(read.table(list.files(pattern = "_hisat2_FC$")[1], sep = "\t", skip = 1, header = TRUE))

# Read one file and return its "sorted.bam" column, renamed to the file prefix
read_sorted_bam_column <- function(file_path) {
  prefix <- sub("_hisat2_FC$", "", basename(file_path))  # everything before "_hisat2_FC"
  
  df <- read.table(
    file_path,
    sep = "\t",
    skip = 1,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # pick the column that ends with ".sorted.bam"
  sorted_col <- grep("\\.sorted\\.bam$", names(df), value = TRUE)
  
  out <- df %>% dplyr::select(all_of(sorted_col))
  names(out) <- prefix
  out
}

# List files ending exactly in "_hisat2_FC"
files <- list.files(pattern = "_hisat2_FC$", full.names = TRUE)

# Read Geneid from the FIRST file only
geneid <- read.table(
  files[1],
  sep = "\t",
  skip = 1,
  header = TRUE,
  stringsAsFactors = FALSE
)$Geneid

# Combine columns
final_data <- map_dfc(files, read_sorted_bam_column)

# Set row names
rownames(final_data) <- geneid

# Manually simplify each column name 
colnames(final_data)[colnames(final_data) == "PAR01-001-veh"] <- "PV1"
colnames(final_data)[colnames(final_data) == "PAR01-002-veh"] <- "PV2"
colnames(final_data)[colnames(final_data) == "PAR01-003-veh"] <- "PV3"

colnames(final_data)[colnames(final_data) == "TR01-001-veh"] <- "TV1"
colnames(final_data)[colnames(final_data) == "TR01-002-veh"] <- "TV2"
colnames(final_data)[colnames(final_data) == "TR01-003-veh"] <- "TV3"

colnames(final_data)[colnames(final_data) == "TR01-001-OHT"] <- "TO1"
colnames(final_data)[colnames(final_data) == "TR01-002-OHT"] <- "TO2"
colnames(final_data)[colnames(final_data) == "TR01-003-OHT"] <- "TO3"

colnames(final_data)[colnames(final_data) == "TR01-001-erd"] <- "TE1"
colnames(final_data)[colnames(final_data) == "TR01-002-erd"] <- "TE2"
colnames(final_data)[colnames(final_data) == "TRO1-003-erd"] <- "TE3"

colnames(final_data)[colnames(final_data) == "TR01-001-OHT-erd"] <- "TOE1"
colnames(final_data)[colnames(final_data) == "TR01-002-OHT-erd"] <- "TOE2"
colnames(final_data)[colnames(final_data) == "TR01-003-OHT-erd"] <- "TOE3"

write.table(
  final_data,
  file = "/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/T47D_counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)