library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(fgsea)
library(ggrepel)
library(ggplot2)
library(scales)
library(writexl)
library(corto)
library(DESeq2)
library(ComplexHeatmap)
library(randomForest)
library(cluster)
library(RColorBrewer)
library(colorRamp2)
library(enrichR)
library(msigdbr)
library(dplyr)
library(tidyr)
library(grid)
library(ggtext)
library(showtext)

########### Set seed ##############

set.seed(42)


# Load a Google font that supports Greek
font_add_google("Roboto", "roboto")
showtext_auto()

## Define figure directory ##

data.dir = "/Volumes/T9/Github_T47D_P_TR/T47D_dds_vds/"

final_plotdir = "/Volumes/T9/Github_T47D_P_TR/Final figures/"

### Read in normalized data (...tr contains only TAM-R)

dds_cell_line <- readRDS(paste0(data.dir, "T47D_dds_cell_line.RDS"))

coldata_cell_line = as.data.frame(colData(dds_cell_line))


#### Generate gene annotation (Ensemble ID, gene symbol, biotype)

gene.annot <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(dds_cell_line), keytype = "GENEID", columns = c("SYMBOL","GENEID","GENEBIOTYPE"))

protein.coding.genes = as.vector(gene.annot$GENEID[which(gene.annot$GENEBIOTYPE == "protein_coding")])

##################################################################
##### Retrieve and list Hallmark, Kegg, and Reactome gene sets ###
##################################################################

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

IFN.common = intersect(hallmark.IFNA.genes, hallmark.IFNG.genes)
IFNA.unique <- setdiff(hallmark.IFNA.genes, IFN.common)
IFNG.unique <- setdiff(hallmark.IFNG.genes, IFN.common)

###################################################################
#### 2 Result table (DGEA) #####################################
###################################################################

contrast.TR.PAR = c("cell_line","TR01","PAR")

  dgea.df = as.data.frame(results(dds_cell_line[rownames(dds_cell_line) %in% protein.coding.genes,], contrast.TR.PAR))
  
  # Add gene symbols to result tables 
  
  dgea.df$ID = rownames(dgea.df)
  dgea.df = merge(dgea.df, gene.annot,
                         by.x="ID", by.y="GENEID", 
                         all.x=TRUE, all.y=FALSE  )

###################################################################
#### 3 Result table (GSEA) #####################################
###################################################################

### Generate ranked gene list for each comparison #####
  
  l2fc <- dgea.df$log2FoldChange
  names(l2fc) <- dgea.df$ID
  l2fc <- na.omit(l2fc)
  ranked_genes = sort(l2fc, decreasing = TRUE) # sort the list in decreasing order

# Run fgsea enrichment analysis for HALLMARK genes (remember to double check comparison)

  pathways = HALLMARK_gene_sets_fgsea
  
hallmark_df <- fgseaMultilevel(
    pathways = pathways,
    stats = ranked_genes,
    eps = 0,
    scoreType = "std")


#######################################
########### FIGURES 2025 paper ########
#######################################

# Clean up and convert all pathway names to uppercase
hallmark_df$pathway <- gsub("_", " ", gsub("HALLMARK_", "", hallmark_df$pathway))
hallmark_df$pathway <- toupper(hallmark_df$pathway)

# Filter the significant pathways once and store them
df_sig <- hallmark_df[hallmark_df$padj < 0.05, ]
df_sig <- df_sig[order(df_sig$padj)[1:10], ]

# Boldface specific interferon pathways (matching in uppercase)
df_sig$fontface <- ifelse(df_sig$pathway %in% c("INTERFERON ALPHA RESPONSE", "INTERFERON GAMMA RESPONSE"), "bold", "plain")

# Create label column and apply line breaks in uppercase
df_sig$label <- df_sig$pathway
df_sig$label[df_sig$pathway == "EPITHELIAL MESENCHYMAL TRANSITION"] <- "EPITHELIAL MESENCHYMAL\nTRANSITION"
df_sig$label[df_sig$pathway == "INTERFERON ALPHA RESPONSE"] <- "IFN-\u03B1 RESPONSE"
df_sig$label[df_sig$pathway == "INTERFERON GAMMA RESPONSE"] <- "IFN-\u03B3 RESPONSE"
df_sig$label[df_sig$pathway == "ESTROGEN RESPONSE LATE"] <- "ESTROGEN RESPONSE\nLATE"
df_sig$label[df_sig$pathway == "ESTROGEN RESPONSE EARLY"] <- "ESTROGEN RESPONSE\nEARLY"
df_sig$label[df_sig$pathway == "UNFOLDED PROTEIN RESPONSE"] <- "UNFOLDED PROTEIN\nRESPONSE"

# Plot
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
  scale_x_continuous(limits = c(-2.5, 2.5))+
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


ggsave(paste0(final_plotdir,"T47D TR-T47D vs P-T47D GSEA barplot.pdf"), plot = gsea_plot, height = 5, width = 5.93, units = "cm")


########################################################################
########## Figure 1D Volcano plot of TR-T47D vs. P-T47D      ########### 
########################################################################

# Diffeexpressed cutoff
vul.pcutoff = 5*10^-2
vul.logcutoff1 = 0
vul.logcutoff2 = 0

### Labelcutoff
label.pcutoff = 5*10^-2
label.logcutoff1 = 0
label.logcutoff2 = 0

GSOI1 = IFNA.unique
GSOI1.title = "IFN alpha response genes"
GSOI2 = IFNG.unique
GSOI2.title = "IFN gamma response genes"
GSOI3 = IFN.common
GSOI3.title = "Intersecting IFN alpha/gamma response genes"

## Get contrast of interest ###
DOI = dgea.df

# add a column of NAs
DOI$diffexpressed <- "NO"
# if log2Foldchange > vul.logcutoff2 and padj < vul.pcutoff, set as "UP" 
DOI$diffexpressed[DOI$log2FoldChange > vul.logcutoff2 & DOI$padj < vul.pcutoff] <- "UP"
# if log2Foldchange < vul.logcutoff1 and padj < vul.pcutoff, set as "DOWN"
DOI$diffexpressed[DOI$log2FoldChange < vul.logcutoff1 & DOI$padj < vul.pcutoff] <- "DOWN"

# if gene is within GSOI, paste GSOI title
DOI$GSOI = NA
DOI$GSOI[DOI$SYMBOL %in% GSOI1] <- paste0(GSOI1.title)
DOI$GSOI[DOI$SYMBOL %in% GSOI2] <- paste0(GSOI2.title)
DOI$GSOI[DOI$SYMBOL %in% GSOI3] <- paste0(GSOI3.title)

# add a column of NAs
DOI$delabel <- NA
# if log2Foldchange > label.logcutoff2 and padj < label.pcutoff, set as "UP" 
DOI$delabel[DOI$log2FoldChange > label.logcutoff2 & DOI$padj < label.pcutoff & !is.na(DOI$padj) & !is.na(DOI$GSOI)] <- DOI$SYMBOL[DOI$log2FoldChange > label.logcutoff2 & DOI$padj < label.pcutoff & !is.na(DOI$padj) & !is.na(DOI$GSOI)]
# if log2Foldchange < label.logcutoff1 and padj < label.pcutoff, set as "DOWN"
DOI$delabel[DOI$log2FoldChange < label.logcutoff1 & DOI$padj < label.pcutoff & !is.na(DOI$padj) & !is.na(DOI$GSOI)] <- DOI$SYMBOL[DOI$log2FoldChange < label.logcutoff1 & DOI$padj < label.pcutoff & !is.na(DOI$padj) & !is.na(DOI$GSOI)]
# Alpha differentially expressed genes
DOI$alpha <- 1
DOI$alpha[!is.na(DOI$GSOI) & DOI$diffexpressed != "NO"] <- 2

vol.DOI <- ggplot(data = DOI, aes(x = log2FoldChange, y = as.numeric(-log10(padj)), label = delabel)) +
  ggtitle("TR-T47D vs. P-T47D") +
  ylab("-log10(FDR)") + 
  xlab("log2FoldChange") +
  
  # Background points
  geom_point(data = DOI[is.na(DOI$GSOI), ], size = 0.50, color = "#E0E0E0") +
  
  # Highlighted genes using shape and fixed color
  geom_point(data = DOI[!is.na(DOI$GSOI), ],
             aes(shape = GSOI, alpha = alpha),
             color = "#f28e2b", size = 1) +
  
  scale_alpha(range = c(0.5, 0.9)) +
  geom_text_repel(segment.size = 0.1, cex = 2, max.overlaps = 25) +
  
  # Custom shapes
  scale_shape_manual(
    values = c(
      "IFN alpha response genes" = 16,       # circle
      "Intersecting IFN alpha/gamma response genes" = 17,  # triangle
      "IFN gamma response genes" = 15        # square
    ),
    labels = c(
      expression("IFN-" * alpha * " response genes"),
      expression("IFN-" * gamma * " response genes"),
      expression("Intersecting IFN-" * alpha * "/" * gamma * " response genes")
    )
  ) +
  
  geom_hline(yintercept = -log10(5 * 10^-2), col = "grey39", linetype = 2, linewidth = 0.25) +
  scale_y_continuous(trans = "log1p") +
  scale_x_continuous() +
  theme(
    axis.ticks = element_line(linewidth = 0.25, colour = "black"),
    legend.key = element_rect(fill = NA, color = NA),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),           # vertical spacing between items
    legend.key.height = unit(3, "pt"),          # height of each key (symbol row)
    legend.box.spacing = unit(0, "pt"),
    plot.title = element_blank(),
    legend.position = "top",
    plot.margin = ggplot2::margin(c(0.4,0.2,0.2,0.2), unit = "cm"),
    panel.background = element_blank(),
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6),
    title = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(1, "mm"),
    axis.line = element_line(linewidth = 0.25, colour = "black")
  ) +
  guides(shape = guide_legend(ncol = 1), alpha = "none")

vol.DOI

ggsave(
  filename = paste0(final_plotdir, "T47D TR-T47D vs P-T47D volcano plot.pdf"),
  plot = vol.DOI,
  width = 59.26,
  height = 55,
  units = "mm"
)


