library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(fgsea)
library(ggrepel)
library(ggplot2)
library(ggsignif)  
library(dplyr)
library(rlang)
library(showtext)

# Load a Google font that supports Greek
font_add_google("Roboto", "roboto")
showtext_auto()

## Define figure directory ##

data.dir = "/Volumes/T9/Github_projects/FUTURE_FGFR1_STING/Github_T47D_P_TR/data/"

final_plotdir = "/Volumes/T9/Github_T47D_P_TR_figures/Final figures/"

supp_plotdir = "/Volumes/T9/Github_T47D_P_TR_figures/Final supplementary figures/"

### Read in normalized data 

dds_treatment <- readRDS(paste0(data.dir, "T47D_dds_treatment.RDS"))

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

# ERD vs. CTR
make_fixed_element_gsea_plot(
  hallmark_df = gsea.list$HALLMARK$hallmark.ERD.CTR,
  plot_title = "erdafitinib vs. control",
  file_path = paste0(final_plotdir, "Fig 2C (TR-T47D erdafitinib vs. control GSEA barplot).pdf")
)

# TE vs. ERD
make_fixed_element_gsea_plot(
  hallmark_df = gsea.list$HALLMARK$hallmark.TE.ERD,
  plot_title = "tamoxifen + erdafitinib vs. erdafitinib",
  file_path = paste0(final_plotdir, "Fig 2E (TR-T47D tamoxifen + erdafitinib vs. erdafitinib GSEA barplot).pdf")
)


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
  output_file = paste0(final_plotdir, "Fig 2D (TR-T47D erdafitinib vs. control volcano plot).pdf")
)

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.TAM.CTR,
  geneset = hallmark.IFN.genes,
  plot_title = "4-OHT vs. control",
  output_file = paste0(supp_plotdir, "Fig S2E (TR-T47D 4-OHT vs. control volcano plot).pdf")
)

plot_volcano_highlight_geneset(
  dgea_df = dgea.list$dgea.TE.ERD,
  geneset = hallmark.IFN.genes,
  plot_title = "4-OHT + erdafitinib vs. erdafitinib",
  output_file = paste0(final_plotdir, "Fig 2F (TR-T47D 4-OHT + erdafitinib vs. erdafitinib volcano plot).pdf")
)

