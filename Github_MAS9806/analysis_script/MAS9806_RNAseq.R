library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(fgsea)
library(ggrepel)
library(ggplot2)
library(msigdbr)
library(dplyr)
library(showtext)

########### Set seed ##############

set.seed(42)

# Load a Google font that supports Greek
font_add_google("Roboto", "roboto")
showtext_auto()

## Define figure directory ##

final.data.dir = "/Volumes/T9/Github_projects/Future_FGFR1_STING/Github_MAS9806/data/"

final_plotdir = "/Volumes/T9/Github_MAS9806_figures/Final figures/"

supp_plotdir = "/Volumes/T9/Github_MAS9806_figures/Final supplementary figures/"

supp_tabledir = "/Volumes/T9/Github_MAS9806_figures/Final supplementary tables/"

### Read in normalized data 

dds <- readRDS(paste0(final.data.dir, "MAS98.06_dds.RDS")) ## Script for dds/vsd files are in Final data/MAS9806_RNAseq/

#### Generate gene annotation (Ensemble ID, gene symbol, biotype)

gene.annot <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(dds), keytype = "GENEID", columns = c("SYMBOL","GENEID","GENEBIOTYPE"))

protein.coding.genes = as.vector(gene.annot$GENEID[which(gene.annot$GENEBIOTYPE == "protein_coding")])


##################################################################
##### Retrieve and list Hallmark, Kegg, and Reactome gene sets ###
##################################################################

### Retrieve HALLMARK gene sets for Homo sapiens from MsigDB (database)
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)

HALLMARK_gene_sets_fgsea = split(x = HALLMARK_gene_sets$ensembl_gene, f = HALLMARK_gene_sets$gs_name)

#########################################################################################################
#### Generate character vectors of gene symbols/ensemble IDs within indivudual hallmark gene sets #######
#########################################################################################################

# Hallmark Interferon response alpha/gamma
hallmark.IFNA.genes = HALLMARK_gene_sets$gene_symbol[HALLMARK_gene_sets$gs_name =="HALLMARK_INTERFERON_ALPHA_RESPONSE"]
hallmark.IFNG.genes = HALLMARK_gene_sets$gene_symbol[HALLMARK_gene_sets$gs_name =="HALLMARK_INTERFERON_GAMMA_RESPONSE"]
hallmark.IFN.genes = union(hallmark.IFNA.genes,hallmark.IFNG.genes)

###################################################################
#### 2 Result tables (DGEA) #####################################
###################################################################

dgea.list <- list()
cont.list <- list()

# Define contrasts and assign names to cont.list
cont.list[["TAM.CTR"]] <- c("treatment", "TAM", "CTR")
cont.list[["ERD.CTR"]] <- c("treatment", "ERD", "CTR")
cont.list[["TE.CTR"]]  <- c("treatment", "TE",  "CTR")
cont.list[["TE.TAM"]]  <- c("treatment", "TE",  "TAM")
cont.list[["TE.ERD"]]  <- c("treatment", "TE",  "ERD")
cont.list[["TAM.ERD"]] <- c("treatment", "TAM", "ERD")

# Iterate through contrasts and compute results
for (name in names(cont.list)) {
  
  # Filter based on protein-coding genes
  res <- results(dds[rownames(dds) %in% protein.coding.genes, ], cont.list[[name]])
  
  # Convert to data frame and add gene symbols
  res <- as.data.frame(res)
  res$ID <- rownames(res)
  res <- merge(res, gene.annot, by.x = "ID", by.y = "GENEID", all.x = TRUE, all.y = FALSE)
  
  # Store in dgea.list with proper naming
  dgea.list[[paste0("dgea.", name)]] <- res
}

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

HALLMARK_gsea_list = list()


# Run fgsea enrichment analysis for HALLMARK genes (remember to double check comparison)
pathways = HALLMARK_gene_sets_fgsea

for(i in 1:length(dgea.list)){
  HALLMARK_gsea_list[[i]] <- fgseaMultilevel(
    pathways = pathways,
    stats = ranked_gene_list[[i]],
    eps = 0,
    scoreType = "std"
  )
}

names(HALLMARK_gsea_list) <- gsub("dgea\\.","hallmark.", names(dgea.list))

##########################################################################################
#################################### 2025 IFN paper figures     ###########################
##########################################################################################

make_fixed_element_gsea_plot <- function(gsea_df, plot_title, file_path) {
  
  
  # Clean and format pathway names
  gsea_df$pathway <- gsub("_", " ", gsub("HALLMARK_", "", gsea_df$pathway))
  gsea_df$pathway <- toupper(gsea_df$pathway)
  
  # Filter for significant pathways (padj < 0.05)
  df_sig <- gsea_df %>% filter(padj < 0.05)
  if (nrow(df_sig) == 0) {
    warning("No significant pathways (padj < 0.05), skipping plot.")
    return(NULL)
  }
  
  # Order by padj and take top 5 from +NES and top 5 from -NES
  df_sig_neg <- df_sig[df_sig$NES < 0, ]
  df_sig_pos <- df_sig[df_sig$NES > 0, ]
  
  top_neg <- df_sig_neg[order(df_sig_neg$padj)[1:min(5, nrow(df_sig_neg))], ]
  top_pos <- df_sig_pos[order(df_sig_pos$padj)[1:min(5, nrow(df_sig_pos))], ]
  
  df_sig <- rbind(top_neg, top_pos)
  
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
  df_sig$label[df_sig$pathway == "OXIDATIVE PHOSPHORYLATION"] <- "OXIDATIVE\nPHOSPHORYLATION"
  
  
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
    scale_x_continuous(limits = c(-3.5, 3.5))+
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
      plot.margin = unit(c(0.4, 0.2, 0.1, 0.2), "cm")
    )
  
  # Save the figure
  ggsave(file_path, plot = gsea_plot, height = 53, width = 59.26, units = "mm")
  
  return(gsea_plot)
}

######### GSEA barplots ##############

# TAM vs. CTR
make_fixed_element_gsea_plot(
  gsea_df = HALLMARK_gsea_list$hallmark.TAM.CTR,
  plot_title = "tamoxifen vs. control",
  file_path = paste0(final_plotdir, "MAS9806 tamoxifen vs. control GSEA barplot.pdf")
)

# ERD vs. CTR
make_fixed_element_gsea_plot(
  gsea_df = HALLMARK_gsea_list$hallmark.ERD.CTR,
  plot_title = "erdafitinib vs. control",
  file_path = paste0(final_plotdir, "MAS9806 erdafitinib vs. control GSEA barplot.pdf")
)

# TE vs. TAM
make_fixed_element_gsea_plot(
  gsea_df = HALLMARK_gsea_list$hallmark.TE.TAM,
  plot_title = "tamoxifen + erdafitinib vs. tamoxifen",
  file_path = paste0(final_plotdir, "MAS9806 tamoxifen + erdafitinib vs. tamoxifen GSEA barplot.pdf")
)

##########################################################################################
#################################### Volcano plots             ###########################
##########################################################################################

## Volcano plot function
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
    highlight_color = "#F28E2B",
    geneset_label = "Highlighted genes"  # NEW ARGUMENT
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
    
    # Highlighted genes (using dynamic legend label)
    geom_point(
      data = df[df$SYMBOL %in% geneset, ],
      aes(color = factor(geneset_label), alpha = alpha),
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
    scale_color_manual(values = setNames(highlight_color, geneset_label)) +
    
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

## Volcano plots of various comparisons highlighting IFN response genes

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.ERD.CTR,
  geneset = hallmark.IFN.genes,
  geneset_label = "IFN response genes",
  output_file = paste0(final_plotdir, "MAS98.06 erdafitinib vs. control volcano plot.pdf")
)

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.TAM.CTR,
  geneset = hallmark.IFN.genes,
  geneset_label = "IFN response genes",
  output_file = paste0(final_plotdir, "MAS98.06 tamoxifen vs. control volcano plot test.pdf")
)

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.TE.TAM,
  geneset = hallmark.IFN.genes,
  geneset_label = "IFN response genes",
  output_file = paste0(final_plotdir, "MAS98.06 tamoxifen + erdafitinib vs. tamoxifen volcano plot.pdf")
)

##########################################
## Volcano plot of FGFR1 + STING #########
##########################################

plot_volcano_custom_highlight <- function(
    dgea_df,
    highlight_genes = NULL,  # Named vector: names are gene symbols, values are custom labels
    highlight_colors = NULL, # Named vector: names are gene symbols, values are custom colors
    plot_title = "Volcano Plot",
    output_file = NULL,
    vul_pcutoff = 0.05,
    vul_logcutoff1 = 0,
    vul_logcutoff2 = 0
) {
  df <- dgea_df
  
  # Differential expression classification
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange > vul_logcutoff2 & df$padj < vul_pcutoff] <- "UP"
  df$diffexpressed[df$log2FoldChange < vul_logcutoff1 & df$padj < vul_pcutoff] <- "DOWN"
  
  # Label column
  df$label <- NA
  if (!is.null(highlight_genes)) {
    df$label[df$SYMBOL %in% names(highlight_genes)] <- highlight_genes[df$SYMBOL[df$SYMBOL %in% names(highlight_genes)]]
  }
  
  # Alpha and color columns
  df$highlight <- ifelse(df$SYMBOL %in% names(highlight_genes), "highlight", "background")
  df$alpha <- ifelse(df$highlight == "highlight", 2, 1)
  df$color <- ifelse(df$SYMBOL %in% names(highlight_colors), highlight_colors[df$SYMBOL], "#E0E0E0")
  
  # Plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), label = label)) +
    ggtitle(plot_title) +
    ylab("-log10(FDR)") +
    xlab("log2FoldChange") +
    
    # Background points
    geom_point(data = df[df$highlight == "background", ], size = 0.5, color = "#E0E0E0") +
    
    # Highlighted points with color
    geom_point(data = df[df$highlight == "highlight", ],
               aes(color = SYMBOL, alpha = alpha), size = 1) +
    
    # Labels
    geom_text_repel(segment.size = 0.1, cex = 2, max.overlaps = 25) +
    
    # Significance line
    geom_hline(yintercept = -log10(vul_pcutoff), col = "grey39", linetype = 3, linewidth = 0.5) +
    
    # Scales
    scale_alpha(range = c(0.5, 0.9)) +
    scale_y_continuous(trans = "log1p") +
    scale_color_manual(values = highlight_colors) +
    
    # Theme
    theme(
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      legend.position = "none",
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
    guides(alpha = "none")
  
  # Print and optionally save
  print(p)
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 59.26, height = 55, units = "mm")
  }
  
  return(p)
}

############# Volcano plot for STING + FGFR1 in tamoxifen vs. control 

# Define your gene highlights and colors
highlight_genes <- c(TMEM173 = "STING1", FGFR1 = "FGFR1")
highlight_colors <- c(TMEM173 = "#D37295", FGFR1 = "#59A14F")

# Run the plot function
plot_volcano_custom_highlight(
  dgea_df = dgea.list$dgea.TAM.CTR,
  highlight_genes = highlight_genes,
  highlight_colors = highlight_colors,
  output_file = paste0(supp_plotdir, "MAS98.06 tamoxifen vs. control volcano plot_FGFR1_STING.pdf")
)

############################################################
###### Supplementary tables for all DGEA comparisons #######
############################################################

# Load required package
library(openxlsx)

# Define name mapping
sheet_name_map <- c(
  "dgea.TAM.CTR" = "tamoxifen vs. Vehicle",
  "dgea.ERD.CTR" = "Erdafitinib vs. Vehicle",
  "dgea.TE.CTR"  = "tamoxifen+Erd vs. Vehicle",
  "dgea.TE.TAM"  = "tamoxifen+Erd vs. tamoxifen",
  "dgea.TE.ERD"  = "tamoxifen+Erd vs. Erdafitinib",
  "dgea.TAM.ERD" = "tamoxifen vs. Erdafitinib"
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
output_path <- file.path(supp_tabledir, "MAS98.06_DGEA_results_all_comparisons.xlsx")

# Save the workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)

# Confirmation
cat("Workbook saved to:", output_path, "\n")

#####################################################################
###### Supplementary tables for all Hallmark GSEA comparisons #######
#####################################################################

# Define sheet name mapping
gsea_sheet_map <- c(
  "hallmark.TAM.CTR" = "tamoxifen vs. Vehicle",
  "hallmark.ERD.CTR" = "Erdafitinib vs. Vehicle",
  "hallmark.TE.CTR"  = "tamoxifen+Erd vs. Vehicle",
  "hallmark.TE.TAM"  = "tamoxifen+Erd vs. tamoxifen",
  "hallmark.TE.ERD"  = "tamoxifen+Erd vs. Erdafitinib",
  "hallmark.TAM.ERD" = "tamoxifen vs. Erdafitinib"
)

# Create a new workbook
wb <- createWorkbook()

# Loop through GSEA results and write to workbook
for (name in names(HALLMARK_gsea_list)) {
  sheet_label <- gsea_sheet_map[[name]]
  gsea_table <- HALLMARK_gsea_list[[name]] %>%
    arrange(padj)  # sort by adjusted p-value
  addWorksheet(wb, sheetName = sheet_label)
  writeData(wb, sheet = sheet_label, x = gsea_table)
}

# Define output path
output_path <- file.path(supp_tabledir, "MAS98.06_HALLMARK_GSEA_results.xlsx")

# Save the workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)

# Confirmation
cat("Workbook saved to:", output_path, "\n")


