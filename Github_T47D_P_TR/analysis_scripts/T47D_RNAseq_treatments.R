library(stringr)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(fgsea)
library(ggrepel)
library(ggplot2)
library(GSVA)
library(ggsignif)  
library(dplyr)
library(rlang)
library(showtext)

# Load a Google font that supports Greek
font_add_google("Roboto", "roboto")
showtext_auto()

### Load color schemes ####

#colour.df <- read.table("/Users/torbjo2503/Desktop/PhD/2024 IFN paper/colour schemes.txt", header = TRUE, sep = " ")

## Define figure directory ##

data.dir = "/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/"

final_plotdir = "/Volumes/T9/Github_T47D_P_TR/Final figures/"

supp_plotdir = "/Volumes/T9/Github_T47D_P_TR/Final supplementary figures/"

supp_tabledir = "/Volumes/T9/Github_T47D_P_TR/Final supplementary tables/"

### Read in normalized data (...tr contains only TAM-R treatments)

dds_treatment <- readRDS(paste0(data.dir, "T47D_dds_treatment.RDS"))
vsd_treatment <- readRDS(paste0(data.dir, "T47D_vsd_treatment.RDS"))

coldata_treatments = as.data.frame(colData(dds_treatment))

########### Set seed ##############

set.seed(42)

#### Generate gene annotation (Ensemble ID, gene symbol, biotype)

gene.annot <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(dds_treatment), keytype = "GENEID", columns = c("SYMBOL","GENEID","GENEBIOTYPE"))

protein.coding.genes = as.vector(gene.annot$GENEID[which(gene.annot$GENEBIOTYPE == "protein_coding")])


### Retrieve HALLMARK gene sets for Homo sapiens from MsigDB (database)
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)

HALLMARK_gene_sets_fgsea = split(x = HALLMARK_gene_sets$ensembl_gene, f = HALLMARK_gene_sets$gs_name)

HALLMARK_gene_sets_symbols_fgsea = split(x = HALLMARK_gene_sets$gene_symbol , f = HALLMARK_gene_sets$gs_name)


#########################################################################################################
#### Generate character vectors of gene symbols/ensemble IDs within indivudual hallmark gene sets #######
#########################################################################################################

# Hallmark Interferon response alpha/gamma
hallmark.IFNA.genes = HALLMARK_gene_sets$gene_symbol[HALLMARK_gene_sets$gs_name =="HALLMARK_INTERFERON_ALPHA_RESPONSE"]
hallmark.IFNG.genes = HALLMARK_gene_sets$gene_symbol[HALLMARK_gene_sets$gs_name =="HALLMARK_INTERFERON_GAMMA_RESPONSE"]
hallmark.IFN.genes = union(hallmark.IFNA.genes,hallmark.IFNG.genes)

hallmark.IFNA.ensemble = HALLMARK_gene_sets$ensembl_gene[HALLMARK_gene_sets$gs_name =="HALLMARK_INTERFERON_ALPHA_RESPONSE"]
hallmark.IFNG.ensemble = HALLMARK_gene_sets$ensembl_gene[HALLMARK_gene_sets$gs_name =="HALLMARK_INTERFERON_GAMMA_RESPONSE"]
hallmark.IFN.ensemble = union(hallmark.IFNA.ensemble,hallmark.IFNG.ensemble)

TISG.genes = c("SAMHD1","OAS1","OAS2","OAS3","DDX60","XAF1","IFI27","IFI44","PLSCR1","HERC6","PARP9","RSAD2","IFIT1","STAT1","MX1","DTX3L","B2M")
TISG.ensemble = gene.annot$GENEID[gene.annot$SYMBOL %in% TISG.genes]

###################################################################
#### 2 Result tables (DGEA) #####################################
###################################################################

dgea.list = list()
cont.list = list()

contrast.TAM.CTR = c("treatment","TAM","CTR")
cont.list[[1]] = contrast.TAM.CTR
contrast.ERD.CTR = c("treatment","ERD","CTR")
cont.list[[2]] = contrast.ERD.CTR
contrast.TE.CTR = c("treatment","TE","CTR")
cont.list[[3]] = contrast.TE.CTR
contrast.TE.TAM = c("treatment","TE","TAM")
cont.list[[4]] = contrast.TE.TAM
contrast.TE.ERD = c("treatment","TE","ERD")
cont.list[[5]] = contrast.TE.ERD
contrast.TAM.ERD = c("treatment","TAM","ERD")
cont.list[[6]] = contrast.TAM.ERD

for(i in 1:6){
  
  # Filter out genes based on protein-coding biotype
  dgea.list[[i]] = as.data.frame(results(dds_treatment[rownames(dds_treatment) %in% protein.coding.genes,], cont.list[[i]]))
  
  # Add gene symbols to result tables 
  
  dgea.list[[i]]$ID = rownames(dgea.list[[i]])
  dgea.list[[i]] = merge(dgea.list[[i]], gene.annot,
                         by.x="ID", by.y="GENEID", 
                         all.x=TRUE, all.y=FALSE  )
  
}

names(dgea.list) <- c("dgea.TAM.CTR","dgea.ERD.CTR","dgea.TE.CTR","dgea.TE.TAM","dgea.TE.ERD", "dgea.TAM.ERD")


###################################################################
#### 3 Result tables (GSEA) #####################################
###################################################################

### Generate ranked gene list for each comparison #####

ranked_gene_list = list()

for(i in 1:length(dgea.list)){
  l2fc <- dgea.list[[i]]$log2FoldChange
  names(l2fc) <- dgea.list[[i]]$ID
  l2fc <- na.omit(l2fc)
  ranked_gene_list[[i]] = sort(l2fc, decreasing = TRUE) # sort the list in decreasing order
  names(ranked_gene_list)[i] = names(dgea.list)[i]
}

### Generate list of GSEA result tables

gsea.hallmark = list()
gsea.kegg = list()
gsea.reactome = list()
gsea.gobp = list()

gsea.list = list(gsea.hallmark, gsea.kegg, gsea.reactome, gsea.gobp)

names(gsea.list) = c("HALLMARK","KEGG","REACTOME", "GOBP")


# Run fgsea enrichment analysis for HALLMARK genes (remember to double check comparison)
pathways = HALLMARK_gene_sets_fgsea

for(i in 1:length(dgea.list)){
  gsea.list[["HALLMARK"]][[i]] <- fgseaMultilevel(
    pathways = pathways,
    stats = ranked_gene_list[[i]],
    eps = 0,
    scoreType = "std"
  )
}

names(gsea.list[["HALLMARK"]]) <- gsub("dgea\\.","hallmark.", names(dgea.list))


#######################################
########### FIGURES  ########
#######################################

### GSEA barplot function

make_fixed_element_gsea_plot <- function(hallmark_df, plot_title, file_path) {

  
  # Clean and format pathway names
  hallmark_df$pathway <- gsub("_", " ", gsub("HALLMARK_", "", hallmark_df$pathway))
  hallmark_df$pathway <- toupper(hallmark_df$pathway)
  
  # Filter for significant pathways (padj < 0.05)
  df_sig <- hallmark_df %>% filter(padj < 0.05)
  if (nrow(df_sig) == 0) {
    warning("No significant pathways (padj < 0.05), skipping plot.")
    return(NULL)
  }
  
  # Order by padj and take top 10
  df_sig <- df_sig[order(df_sig$padj)[1:min(10, nrow(df_sig))], ]
  
  # Set bold fontface for interferon responses
  df_sig$fontface <- ifelse(df_sig$pathway %in% c("INTERFERON ALPHA RESPONSE", "INTERFERON GAMMA RESPONSE"),
                            "bold", "plain")
  
  # Apply manual line breaks to known long names
  df_sig$label <- df_sig$pathway
  df_sig$label[df_sig$pathway == "EPITHELIAL MESENCHYMAL TRANSITION"] <- "EPITHELIAL MESENCHYMAL\nTRANSITION"
  df_sig$label[df_sig$pathway == "INTERFERON ALPHA RESPONSE"] <- "IFN-\u03B1 RESPONSE"
  df_sig$label[df_sig$pathway == "INTERFERON GAMMA RESPONSE"] <- "IFN-\u03B3 RESPONSE"
  df_sig$label[df_sig$pathway == "ESTROGEN RESPONSE LATE"] <- "ESTROGEN RESPONSE\nLATE"
  df_sig$label[df_sig$pathway == "ESTROGEN RESPONSE EARLY"] <- "ESTROGEN RESPONSE\nEARLY"
  df_sig$label[df_sig$pathway == "UNFOLDED PROTEIN RESPONSE"] <- "UNFOLDED PROTEIN\nRESPONSE"
  df_sig$label[df_sig$pathway == "ALLOGRAFT REJECTION"] <- "ALLOGRAFT\nREJECTION"
  df_sig$label[df_sig$pathway == "CHOLESTEROL HOMEOSTASIS"] <- "CHOLESTEROL\nHOMEOSTASIS"
  df_sig$label[df_sig$pathway == "INFLAMMATORY RESPONSE"] <- "INFLAMMATORY\nRESPONSE"
  
  
  
  # Create the plot
  gsea_plot <- ggplot(df_sig, aes(x = NES, y = reorder(pathway, NES), 
                                  fill = ifelse(NES > 0, "up", "down"))) +
    geom_bar(stat = "identity", width = 0.9) +
    ylab(NULL)+
    geom_text(data = subset(df_sig, NES > 0),
              aes(label = label, fontface = fontface),
              x = -0.05, hjust = 1, nudge_x = 0.1,
              color = "black", size = 5 / .pt,
              lineheight = 0.8) +
    
    geom_text(data = subset(df_sig, NES < 0),
              aes(label = label, fontface = fontface),
              x = 0.05, hjust = 0, nudge_x = -0.1,
              color = "black", size = 5 / .pt,
              lineheight = 0.8) +
    
    scale_fill_manual(values = c("up" = "#e15759", "down" = "#4e79a7")) +
    labs(x = "Normalized Enrichment Score (NES)", y = "", title = NULL) +
    scale_x_continuous(limits = c(-2.5, 3.5))+
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3),
      axis.text.x = element_text(color = "black", size = 6),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0.4, 0.2, 0.2, 0.2), "cm")
    )
  
  # Save the figure
  ggsave(file_path, plot = gsea_plot, height = 50, width = 59.26, units = "mm")
  
  return(gsea_plot)
}


######### GSEA barplots ##############

# TAM vs. CTR
make_fixed_element_gsea_plot(
  hallmark_df = gsea.list$HALLMARK$hallmark.TAM.CTR,
  plot_title = "tamoxifen vs. control",
  file_path = paste0(supp_plotdir, "TR-T47D tamoxifen vs. control GSEA barplot.pdf")
)

# ERD vs. CTR
make_fixed_element_gsea_plot(
  hallmark_df = gsea.list$HALLMARK$hallmark.ERD.CTR,
  plot_title = "erdafitinib vs. control",
  file_path = paste0(final_plotdir, "TR-T47D erdafitinib vs. control GSEA barplot.pdf")
)

# TE vs. ERD
make_fixed_element_gsea_plot(
  hallmark_df = gsea.list$HALLMARK$hallmark.TE.ERD,
  plot_title = "tamoxifen + erdafitinib vs. erdafitinib",
  file_path = paste0(final_plotdir, "TR-T47D tamoxifen + erdafitinib vs. erdafitinib GSEA barplot.pdf")
)

# TE vs. TAM
#make_fixed_element_gsea_plot(
  #hallmark_df = gsea.list$HALLMARK$hallmark.TE.TAM,
 # plot_title = "tamoxifen + erdafitinib vs. tamoxifen",
 # file_path = paste0(supp_plotdir, "TR-T47D tamoxifen + erdafitinib vs. tamoxifen GSEA barplot.pdf")
#)


############################################################################
################### Volcano plots#####################################
############################################################################


plot_volcano_highlight_geneset <- function(
    dgea_df,
    geneset,
    plot_title = "Volcano Plot",
    output_file = NULL,
    vul_pcutoff = 0.05,
    vul_logcutoff1 = 0,
    vul_logcutoff2 = 0,
    label_pcutoff = 0.05,
    label_logcutoff1 = 0,
    label_logcutoff2 = 0,
    highlight_color = "#F28E2B"
) {

  df <- dgea_df
  
  # Differential expression classification
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange > vul_logcutoff2 & df$padj < vul_pcutoff] <- "UP"
  df$diffexpressed[df$log2FoldChange < vul_logcutoff1 & df$padj < vul_pcutoff] <- "DOWN"
  
  # Label column: only label genes in geneset if they pass significance + FC cutoff
  df$label <- NA
  sig_mask <- df$SYMBOL %in% geneset & df$padj < label_pcutoff & !is.na(df$padj)
  df$label[sig_mask & df$log2FoldChange > label_logcutoff2] <- df$SYMBOL[sig_mask & df$log2FoldChange > label_logcutoff2]
  df$label[sig_mask & df$log2FoldChange < label_logcutoff1] <- df$SYMBOL[sig_mask & df$log2FoldChange < label_logcutoff1]
  
  # Alpha channel for highlighted genes
  df$alpha <- ifelse(df$SYMBOL %in% geneset & df$diffexpressed != "NO", 2, 1)
  
  # Plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), label = label)) +
    ggtitle(plot_title) +
    ylab("-log10(FDR)") +
    xlab("log2FoldChange") +
    
    # Background points
    geom_point(data = df[!df$SYMBOL %in% geneset, ], size = 0.5, color = "#E0E0E0") +
    
    # Highlighted genes (uniform color/shape)
    geom_point(
      data = df[df$SYMBOL %in% geneset, ],
      aes(color = factor("Interferon response genes"), alpha = alpha),
      size = 1
    ) +
    
    # Labels
    geom_text_repel(segment.size = 0.1, cex = 2, max.overlaps = 25) +
    
    # Significance line
    geom_hline(yintercept = -log10(vul_pcutoff), col = "grey39", linetype = 3, linewidth = 0.5) +
    
    # Axis scales
    scale_alpha(range = c(0.5, 0.9)) +
    scale_y_continuous(trans = "log1p") +
    scale_x_continuous() +
    scale_color_manual(values = c("Interferon response genes" = highlight_color), labels = c("Interferon response genes" = "IFN response genes")) +
  
    # Theme
    theme(
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 6),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.y = unit(0, "pt"), 
      legend.key.height = unit(3, "pt"),  
      legend.box.spacing = unit(0, "pt"),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 6),
      plot.title = element_blank(),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      plot.margin = margin(0.4, 0.2, 0.2, 0.2, "cm")
    ) +
    guides(alpha = "none", color = guide_legend(title = NULL))
  
  
  # Print and optionally save
  print(p)
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 59.26, height = 55, units = "mm")
  }
  
  return(p)
}


## Volcano plots of various comparisons

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.ERD.CTR,
  geneset = hallmark.IFN.genes,
  plot_title = "erdafitinib vs. control",
  output_file = paste0(final_plotdir, "TR-T47D erdafitinib vs. control volcano plot.pdf")
)

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.TAM.CTR,
  geneset = hallmark.IFN.genes,
  plot_title = "4-OHT vs. control",
  output_file = paste0(supp_plotdir, "TR-T47D 4-OHT vs. control volcano plot.pdf")
)

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.TE.ERD,
  geneset = hallmark.IFN.genes,
  plot_title = "4-OHT + erdafitinib vs. erdafitinib",
  output_file = paste0(final_plotdir, "TR-T47D 4-OHT + erdafitinib vs. erdafitinib volcano plot.pdf")
)

#plot_volcano_highlight_geneset(
 # dgea_df = dgea.list$dgea.TE.TAM,
 # geneset = hallmark.IFN.genes,
 # plot_title = "4-OHT + erdafitinib vs. 4-OHT",
#  output_file = paste0(supp_plotdir, "TR-T47D 4-OHT + erdafitinib vs. 4-OHT volcano plot.pdf")
#)




########## ssGSEA analysis #################

#--- 1. Prepare: Utilities ---------------------------------------------------#
p_to_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("NS")
}

#--- 2. Define gene set ------------------------------------------------------#
TISG <- c("OAS1", "SAMHD1", "OAS3", "OAS2", "STAT1", "XAF1", "RSAD2", "DDX60",
          "IFI44", "PARP9", "HERC6", "MX1", "DTX3L", "IFI27", "B2M", "IFIT1", "PLSCR1")


# Step 1: Query genes for all TISG symbols
TISG_ensembl_all <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  keys = TISG,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GENEID", "GENEBIOTYPE")
)

# Step 2: Get additional gene metadata including seqnames
gene_metadata <- genes(
  EnsDb.Hsapiens.v86,
  filter = GeneIdFilter(TISG_ensembl_all$GENEID),
  return.type = "data.frame"
)

# Step 3: Keep only genes on standard chromosomes (e.g., chr1–22, X, Y)
standard_chr_genes <- gene_metadata %>%
  filter(seq_name %in% c(as.character(1:22), "X", "Y"))

# Step 4: Merge to get SYMBOLs and deduplicate
TISG_ensembl_canonical <- TISG_ensembl_all %>%
  semi_join(standard_chr_genes, by = c("GENEID" = "gene_id")) %>%
  distinct(SYMBOL, .keep_all = TRUE)

#--- 3. Run ssGSEA -----------------------------------------------------------#
vst_matrix <- assay(vsd_treatment)
custom_gene_set <- list(TISG = TISG_ensembl_canonical$GENEID)
ssgsea_scores <- gsva(ssgseaParam(vst_matrix, custom_gene_set))

#--- 4. Prepare metadata and plot data ---------------------------------------#
plotdata <- cbind(coldata_treatments, t(ssgsea_scores)[rownames(coldata_treatments), , drop = FALSE])
plotdata$treatment <- factor(plotdata$treatment, levels = c("CTR", "TAM", "ERD", "TE"))

#--- 5. Statistics for annotation --------------------------------------------#
comparisons <- list(c("CTR", "TAM"), c("CTR", "ERD"), c("ERD", "TE"), c("CTR", "TE"))
names(comparisons) <- sapply(comparisons, function(x) paste(x, collapse = "_vs_"))

p_values <- sapply(comparisons, function(x) {
  t.test(plotdata[["TISG"]][plotdata$treatment == x[1]],
         plotdata[["TISG"]][plotdata$treatment == x[2]], var.equal = FALSE)$p.value
})

p_adj_values <- p.adjust(p_values, method = "BH")

signif_annotations <- sapply(p_adj_values, p_to_stars)
y_positions <- c(1.0, 1.1, 1.2, 1.3)  # Matching order of comparisons

#--- 6. Plot aesthetics ------------------------------------------------------#
manual_colors <- c("CTR" = "#A0CBE8", "TAM" = "#F1CE63", "ERD" = "#8CD17D", "TE" = "#E15759")
custom_labels <- c("CTR" = "DMSO", "TAM" = "4-OHT", "ERD" = "erdafitinib"
,"TE" = "4-OHT + 
erdafitinib")

#--- 7. Plot generation ------------------------------------------------------#

plot <- ggplot(plotdata, aes(x = treatment, y = TISG, fill = treatment)) +
  ylim(0,1.5)+
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.2) +
  geom_jitter(width = 0.2, size = 0.3, color = "black") +
  geom_signif(comparisons = comparisons,
              annotations = signif_annotations,
              y_position = y_positions,
              tip_length = 0.01,
              vjust = 0,
              textsize = 6 / 2.835,
              map_signif_level = FALSE,
              size = 0.2) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = manual_colors) +
  labs(x = NULL, y = "ssGSEA Enrichment Score", title = "TR-T47D") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5, color = "black", size = 6),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(size = 6),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.position = "none",
    plot.margin = unit(c(0.4, 0.2, 0.2, 0.4), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

#--- 8. Save ------------------------------------------------------------------#
ggsave(filename = file.path(final_plotdir, "TR-T47D_all_treatments_ssGSEA_boxplot.pdf"),
       plot = plot, width = 5.926, height = 5.5, units = "cm")




############################################################
###### Save csv for STRING analysis ####################
############################################################
# Define output file path
output_file <- file.path(data.dir, "dgea_TE_ERD.csv")

# Save the dataframe as CSV
write.csv(dgea.list$dgea.TE.ERD, file = output_file, row.names = FALSE)


############################################################
###### Supplementary tables for all DGEA comparisons #######
############################################################

# Load required package
library(openxlsx)

# Define name mapping
sheet_name_map <- c(
  "dgea.TAM.CTR" = "4-OHT vs. Vehicle",
  "dgea.ERD.CTR" = "Erdafitinib vs. Vehicle",
  "dgea.TE.CTR"  = "4-OHT+Erd vs. Vehicle",
  "dgea.TE.TAM"  = "4-OHT+Erd vs. 4-OHT",
  "dgea.TE.ERD"  = "4-OHT+Erd vs. Erdafitinib",
  "dgea.TAM.ERD" = "4-OHT vs. Erdafitinib"
)

# Create a new workbook
wb <- createWorkbook()

# Loop and write each DGEA result with sorting
for (name in names(dgea.list)) {
  sheet_label <- sheet_name_map[[name]]
  dgea_table <- dgea.list[[name]] %>%
    arrange(padj)  # sort by adjusted p-value
  addWorksheet(wb, sheetName = sheet_label)
  writeData(wb, sheet = sheet_label, x = dgea_table)
}

# Define output path
output_path <- file.path(supp_tabledir, "TR-T47D_DGEA_results_all_comparisons.xlsx")

# Save the workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)

# Confirmation
cat("Workbook saved to:", output_path, "\n")

#####################################################################
###### Supplementary tables for all Hallmark GSEA comparisons #######
#####################################################################

# Define sheet name mapping
gsea_sheet_map <- c(
  "hallmark.TAM.CTR" = "4-OHT vs. Vehicle",
  "hallmark.ERD.CTR" = "Erdafitinib vs. Vehicle",
  "hallmark.TE.CTR"  = "4-OHT+Erd vs. Vehicle",
  "hallmark.TE.TAM"  = "4-OHT+Erd vs. 4-OHT",
  "hallmark.TE.ERD"  = "4-OHT+Erd vs. Erdafitinib",
  "hallmark.TAM.ERD" = "4-OHT vs. Erdafitinib"
)

# Create a new workbook
wb <- createWorkbook()

# Loop through GSEA results and write to workbook
for (name in names(gsea.list$HALLMARK)) {
  sheet_label <- gsea_sheet_map[[name]]
  gsea_table <- gsea.list$HALLMARK[[name]] %>%
    arrange(padj)  # sort by adjusted p-value
  addWorksheet(wb, sheetName = sheet_label)
  writeData(wb, sheet = sheet_label, x = gsea_table)
}

# Define output path
output_path <- file.path(supp_tabledir, "TR-T47D_HALLMARK_GSEA_results.xlsx")

# Save the workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)

# Confirmation
cat("Workbook saved to:", output_path, "\n")




