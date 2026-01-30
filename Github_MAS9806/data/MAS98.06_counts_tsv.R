library(purrr)

data.dir = ("/Volumes/T9/MAS9806/MAS9806_RNAseq_dds/30_featureCounts/")

# Check that ensemble IDs are equal across all 30 files

setwd(data.dir)
all(sapply(1:30, function(i) {
  all(read.table(sprintf("Future-%d_hisat2_FC", i), sep = "\t", skip = 1, header = TRUE, stringsAsFactors = FALSE)$Geneid ==
        read.table("Future-1_hisat2_FC", sep = "\t", skip = 1, header = TRUE, stringsAsFactors = FALSE)$Geneid)
}))


read_sorted_bam_column <- function(file_number) {
  file_name <- sprintf("Future-%d_hisat2_FC", file_number)
  column_name <- sprintf("Future.%d_hisat2.sorted.bam", file_number)
  
  # Read the file and select the specific column, retaining it as a data frame
  selected_column_df <- read.table(file_name, sep = "\t", skip = 1, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(!!rlang::sym(column_name))
  
  return(selected_column_df)
}

# Apply the function to each file number, combining the results into a single data frame by columns
final_data <- map_dfc(1:30, read_sorted_bam_column)

# Rename all columns to Future-1(1-30)
names(final_data) <- gsub("Future\\.(\\d+)_hisat2\\.sorted\\.bam", "future-\\1", names(final_data))

# Adjusting the function to use read.table
geneid_column <- read.table("Future-1_hisat2_FC", sep = "\t", header = TRUE, skip = 1, stringsAsFactors = FALSE) %>%
  dplyr::pull(Geneid)

# Add the "Geneid" column to the beginning of final_data
final_data <- cbind(Ensembl_ID = geneid_column, final_data)

write.table(final_data,
            file = "/Volumes/T9/MAS9806/MAS9806_RNAseq_dds/mas9806_counts.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
