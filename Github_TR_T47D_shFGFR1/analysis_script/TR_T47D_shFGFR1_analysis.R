library(DESeq2)
library(msigdbr)
library(EnsDb.Hsapiens.v86)
library(fgsea)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)

########### Set seed ##############

set.seed(42)

# Load font with full unicode coverage for greek letters (alpha/beta)

library(showtext)
library(ggtext)

font_add_google("Roboto", "roboto")
showtext_auto()


# Directories

datadir = "/Volumes/T9/Github_projects/Future_FGFR1_STING/Github_TR_T47D_shFGFR1/data/"

final_supp_plotdir = "/Volumes/T9/Github_TR_T47D_shFGFR1_figures/Final supplementary figures/"

# Import dds objects

dds = readRDS(paste0(datadir,"TR_T47D_shFGFR1_dds.rds"))

coldata = as.data.frame(colData(dds))

edb <- EnsDb.Hsapiens.v86

gene_annot <- data.frame(
  gene_id      = rownames(dds),
  gene_name    = mapIds(edb, keys = rownames(dds), column = "GENENAME",
                        keytype = "GENEID", multiVals = "first"),
  gene_biotype = mapIds(edb, keys = rownames(dds), column = "GENEBIOTYPE",
                        keytype = "GENEID", multiVals = "first")
)

gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), ]

pcg_ensemble <- gene_annot$gene_id[!is.na(gene_annot$gene_biotype) & 
                                     gene_annot$gene_biotype == "protein_coding"]
pcg_symbol   <- gene_annot$gene_name[!is.na(gene_annot$gene_biotype) & 
                                       gene_annot$gene_biotype == "protein_coding"]

##################################################################
##### Retrieve Hallmark gene sets ################################
##################################################################

HALLMARK_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets_fgsea <- split(x = HALLMARK_gene_sets$ensembl_gene, f = HALLMARK_gene_sets$gs_name)


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

##################################################################
#### DGEA result table ############################################
##################################################################

dgea_fgfr1_vs_ctr <- as.data.frame(results(dds[rownames(dds) %in% pcg_ensemble,],
                                           contrast = c("treatment", "sh_FGFR1_h2o", "sh_CTR_h2o")))

dgea_fgfr1_vs_ctr$gene_id <- rownames(dgea_fgfr1_vs_ctr)

dgea_fgfr1_vs_ctr <- merge(dgea_fgfr1_vs_ctr, gene_annot[, c("gene_id", "gene_name", "gene_biotype")],
                           by = "gene_id", all.x = TRUE)

##################################################################
#### GSEA on Hallmark ############################################
##################################################################

# Ranked gene list
l2fc <- dgea_fgfr1_vs_ctr$log2FoldChange
names(l2fc) <- dgea_fgfr1_vs_ctr$gene_id
l2fc <- na.omit(l2fc)
ranked_gene_list <- sort(l2fc, decreasing = TRUE)

# Run fgsea
gsea_hallmark <- fgseaMultilevel(
  pathways = HALLMARK_gene_sets_fgsea,
  stats = ranked_gene_list,
  eps = 0,
  scoreType = "std"
)


###############################################################
######################## figures #################
###############################################################

# Prepare GSEA data
df_sig <- as.data.frame(gsea_hallmark)
df_sig$pathway <- gsub("_", " ", gsub("HALLMARK_", "", df_sig$pathway))
df_sig$pathway <- toupper(df_sig$pathway)

df_sig_up   <- df_sig[df_sig$NES > 0 & df_sig$padj < 0.05, ]
df_sig_up   <- df_sig_up[order(df_sig_up$padj)[1:min(5, nrow(df_sig_up))], ]

df_sig_down <- df_sig[df_sig$NES < 0 & df_sig$padj < 0.05, ]
df_sig_down <- df_sig_down[order(df_sig_down$padj)[1:min(5, nrow(df_sig_down))], ]

df_sig <- rbind(df_sig_up, df_sig_down)

df_sig$fontface <- ifelse(df_sig$pathway %in% c("INTERFERON ALPHA RESPONSE",
                                                "INTERFERON GAMMA RESPONSE"),
                          "bold", "plain")

label_map <- c(
  "EPITHELIAL MESENCHYMAL TRANSITION" = "EPITHELIAL MESENCHYMAL\nTRANSITION",
  "INTERFERON ALPHA RESPONSE"         = "IFN-\u03B1 RESPONSE",
  "INTERFERON GAMMA RESPONSE"         = "IFN-\u03B3 RESPONSE",
  "ESTROGEN RESPONSE LATE"            = "ESTROGEN RESPONSE\nLATE",
  "ESTROGEN RESPONSE EARLY"           = "ESTROGEN RESPONSE\nEARLY",
  "UNFOLDED PROTEIN RESPONSE"         = "UNFOLDED PROTEIN\nRESPONSE"
)
df_sig$label <- ifelse(df_sig$pathway %in% names(label_map),
                       label_map[df_sig$pathway],
                       df_sig$pathway)




plot_gsea_hallmark <- function(df_sig,
                               x_limits   = c(-3.2, 3.2),
                               color_up   = "#e15759",
                               color_down = "#4e79a7",
                               save_dir   = NULL,
                               filename   = "gsea_hallmark_barplot") {
  
  p <- ggplot(df_sig, aes(x = NES, y = reorder(pathway, NES),
                          fill = ifelse(NES > 0, "up", "down"))) +
    geom_bar(stat = "identity", width = 0.9) +
    geom_text(data = subset(df_sig, NES > 0),
              aes(label = label, fontface = fontface),
              x = -0.05, hjust = 1, nudge_x = 0.1,
              color = "black", size = 5 / .pt, lineheight = 0.8) +
    geom_text(data = subset(df_sig, NES < 0),
              aes(label = label, fontface = fontface),
              x = 0.05, hjust = 0, nudge_x = -0.1,
              color = "black", size = 5 / .pt, lineheight = 0.8) +
    scale_fill_manual(values = c("up" = color_up, "down" = color_down)) +
    scale_x_continuous(limits = x_limits) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL, title = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      axis.line.x     = element_line(color = "black", linewidth = 0.3),
      axis.ticks.x    = element_line(color = "black", linewidth = 0.3),
      axis.text.x     = element_text(color = "black", size = 6),
      axis.title.x    = element_text(size = 6),
      axis.title.y    = element_blank(),
      axis.line.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.y     = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "none",
      plot.margin     = unit(c(0.4, 0.2, 0.2, 0.2), "cm")
    )
  
  if (!is.null(save_dir)) {
    ggsave(file.path(save_dir, paste0(filename, ".pdf")),
           plot = p, height = 50, width = 59.26, units = "mm")
  }
  
  return(p)
}


gsea_plot <- plot_gsea_hallmark(
  df_sig   = df_sig,
  save_dir = final_supp_plotdir,
  filename = "Fig S2H (TR_T47D_shFGFR1_vs_shCTRL_GSEA_barplot).pdf"
)

#####################################################
################## VOLCANO PLOT #####################
#####################################################

plot_volcano_highlight_geneset <- function(
    dgea_df,
    geneset,
    plot_title      = "Volcano Plot",
    output_file     = NULL,
    vul_pcutoff     = 0.05,
    vul_logcutoff1  = 0,
    vul_logcutoff2  = 0,
    label_pcutoff   = 0.05,
    label_logcutoff1 = 0,
    label_logcutoff2 = 0,
    highlight_color = "#F28E2B"
) {
  
  df <- dgea_df
  
  # Differential expression classification
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange > vul_logcutoff2 & df$padj < vul_pcutoff] <- "UP"
  df$diffexpressed[df$log2FoldChange < vul_logcutoff1 & df$padj < vul_pcutoff] <- "DOWN"
  
  # Label genes in geneset passing significance + FC cutoffs
  df$label <- NA
  sig_mask <- df$gene_name %in% geneset & !is.na(df$padj) & df$padj < label_pcutoff
  df$label[sig_mask & df$log2FoldChange > label_logcutoff2] <- df$gene_name[sig_mask & df$log2FoldChange > label_logcutoff2]
  df$label[sig_mask & df$log2FoldChange < label_logcutoff1] <- df$gene_name[sig_mask & df$log2FoldChange < label_logcutoff1]
  df$alpha <- ifelse(df$gene_name %in% geneset & df$diffexpressed != "NO", 2, 1)
  
  # Plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), label = label)) +
    ylab("-log10(FDR)") +
    xlab("log2FoldChange") +
    geom_point(data = df[!df$gene_name %in% geneset, ], size = 0.5, color = "#E0E0E0") +
    geom_point(data = df[df$gene_name %in% geneset, ],
               aes(color = factor("Interferon response genes"), alpha = alpha), size = 1) +
    geom_text_repel(segment.size = 0.1, cex = 2, max.overlaps = 50) +
    geom_hline(yintercept = -log10(vul_pcutoff), col = "grey39", linetype = 3, linewidth = 0.5) +
    scale_alpha(range = c(0.5, 0.9)) +
    scale_y_continuous(trans = "log1p", breaks = c(0, 2, 5, 10, 20, 30, 50, 100, 200)) +
    scale_x_continuous() +
    scale_color_manual(
      values = c("Interferon response genes" = highlight_color),
      labels = c("Interferon response genes" = "IFN response genes")
    ) +
    theme(
      axis.ticks        = element_line(linewidth = 0.25, color = "black"),
      axis.text         = element_text(size = 6, color = "black"),
      axis.title        = element_text(size = 6),
      axis.line         = element_line(linewidth = 0.25, color = "black"),
      panel.background  = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      plot.margin       = margin(0.4, 0.2, 0.2, 0.2, "cm"),
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.background = element_blank(),
      legend.text       = element_text(size = 6),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.y  = unit(0, "pt"),
      legend.key.height = unit(3, "pt"),
      legend.box.spacing = unit(0, "pt")
    ) +
    guides(alpha = "none", color = guide_legend(title = NULL))
  
  print(p)
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 59.26, height = 55, units = "mm")
  }
  
  return(p)
}

# Call
plot_volcano_highlight_geneset(
  dgea_df     = dgea_fgfr1_vs_ctr,
  geneset     = hallmark.IFN.genes,
  output_file = paste0(final_supp_plotdir, "Fig S2I (TR_T47D_shFGFR1_vs_shCTRL_volcano_plot).pdf")
)


