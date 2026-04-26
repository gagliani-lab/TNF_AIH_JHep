## Figure 2 ##

library(Seurat)
library(dplyr)
library(ggplot2)

AIH.combined <- readRDS(file='./Integration/RDS/3_SCS_SNS_atlas.rds')
hepatocytes.combined <- readRDS(file='./Integration/RDS/14_SNS_hepatocytes_figures.rds')
CD4.combined <- readRDS(file='./Integration/RDS/12_SCS_SNS_CD4_figures.rds')
CD8.combined <- readRDS(file='./Integration/RDS/13_SCS_SNS_CD8_figures.rds')

# Fig 2B -----------------------------------------------------------------------

meta <- AIH.combined@meta.data
meta$UMAP_1 <- AIH.combined@reductions$umap@cell.embeddings[,1]
meta$UMAP_2 <- AIH.combined@reductions$umap@cell.embeddings[,2]
meta <- meta[sample(1:nrow(meta)), ]

ggplot2 <- ggplot(meta, aes(x=UMAP_1,y=UMAP_2,color=protocol)) + 
  geom_point(size=0.01) + 
  scale_color_manual(values=c('#238b45','#cc338b')) +
  theme_classic()
ggsave("SCS_SNS_overlap_umap.jpg", plot = ggplot2, width = 10, height = 10, dpi = 300)

# Fig 2C -----------------------------------------------------------------------

ggplot1 <- DimPlot(immune.combined, label=FALSE, repel = TRUE, 
                   cols = c('CD8' = '#6baed6',
                            'Hepatocytes' = '#006837',
                            'CD4' = '#08519c',
                            'Myeloid' = '#df65b0',
                            'Innate_like_T' = '#1d91c0',
                            'Innate_1' = '#EF6538',
                            'Innate_2' = '#EF6538',
                            'Endothelial' = '#addd8e',
                            'Cholangiocytes' = '#41ab5d',
                            'Stellate/Fibroblasts' = '#807dba', 
                            'B' = '#fd8d3c',
                            'Plasma_cell' = '#cb181d',
                            'Proliferating' = '#c6dbef',
                            'Undefine_1' = '#737373',
                            'Undefine_2' = '#737373',
                            'MT-high' = '#737373' ))
ggsave("all_cell_atlas_umap.jpg", plot = ggplot1, width = 10, height = 10, dpi = 300)

DefaultAssay (immune.combined) <- "RNA"
feature_plots <- FeaturePlot(immune.combined, features = c("PTPRC", "CD3E", "CD4", "CD8A", "ASGR1"),
                             cols = c("lightgrey", "#2171b5"),
                             order = T,
                             slot = "data",
                             min.cutoff = "q5", max.cutoff = "q95",
                             reduction = "umap",
                             pt.size = 0.1,
                             combine= FALSE)
feature_plots

if (is.list(feature_plots)) {
  # Save each plot separately
  for (i in seq_along(feature_plots)) {
    ggsave(paste0("SCS_SNS_atlas_RNA_", feature_plots[[i]]@labels$title, ".jpg"), plot = feature_plots[[i]], width = 10, height = 10, dpi = 300)
  }
} else {
  # If it's a single plot, save it directly
  ggsave(paste0("SCS_SNS_atlas_RNA_", feature_plots@labels$title, ".jpg"), plot = feature_plots, width = 10, height = 10, dpi = 300)
}

DefaultAssay (immune.combined) <- "CITE"
feature_plots_2 <- FeaturePlot(immune.combined, features = c("ADT-CD45", "ADT-TCRAB", "ADT-CD4", "ADT-CD8"),
                               cols = c("lightgrey", "#3f007d"),
                               order = T,
                               slot = "data",
                               min.cutoff = "q5", max.cutoff = "q95",
                               reduction = "umap",
                               pt.size = 0.1,
                               combine= FALSE)
feature_plots_2

if (is.list(feature_plots_2)) {
  # Save each plot separately
  for (i in seq_along(feature_plots_2)) {
    ggsave(paste0("SCS_SNS_atlas_CITE_", feature_plots_2[[i]]@labels$title, ".jpg"), plot = feature_plots_2[[i]], width = 10, height = 10, dpi = 300)
  }
} else {
  # If it's a single plot, save it directly
  ggsave(paste0("SCS_SNS_atlas_CITE_", feature_plots_2@labels$title, ".jpg"), plot = feature_plots_2, width = 10, height = 10, dpi = 300)
}



# Fig 2D -----------------------------------------------------------------------

ggplot1 <- DimPlot(CD4.combined, label= FALSE, pt.size=1.5, group.by="annotation",
                   cols = c('TCM' = '#d99f3e',
                            'TRM1' = '#003c96',
                            'TEM-b' = '#ae98b6',
                            'TREG' = '#807dba',
                            'Naive' = 'lightgrey',
                            'TR1' = '#426e66',
                            'TFH' = '#41ab5d',
                            'TEM-a' = '#c18f6e',
                            'TREG-naive' = '#dadaeb'))
ggsave("SCS_SNS_CD4_umap.jpg", plot = ggplot1, width = 10, height = 10, dpi = 300)

DefaultAssay(CD4.combined) <- "RNA"
genes.heatmap <- c("S1PR1", "CCR7", "SELL", "FAS",
                   "CCR5", "CXCR3",  "ITGAL",
                   "HLA-DRA","HLA-DRB1", "HLA-DQB1",
                   "CCL4","CCL5","CD40LG", "FOS","JUN",
                   "PDCD1","CTLA4","TIGIT","FOXP3","CD8A")
CD4.combined <- ScaleData(CD4.combined,features = genes.heatmap)

avgexp <- as.data.frame(AverageExpression(CD4.combined,features = genes.heatmap,assays = "RNA",slot = "scale.data")$RNA)
avgexp.heatmap <- avgexp[rownames(avgexp) %in% genes.heatmap,] %>% rownames_to_column(var = "gene") %>% gather(-gene,key = "cell.type",value = "avg.expression")
avgexp.heatmap$cell.type = factor(avgexp.heatmap$cell.type, levels = rev( c("Naive", "TCM", "TEM-a", "TEM-b", "TRM1", "TFH", "TR1", "TREG", "TREG-naive")))
avgexp.heatmap$gene = factor(avgexp.heatmap$gene, levels = c("S1PR1", "CCR7", "SELL", "FAS",
                                                             "CXCR3", "CCR5", "ITGAL",
                                                             "HLA-DRA","HLA-DRB1", "HLA-DQB1", 
                                                             "FOS", "JUN", "CD40LG",
                                                             "CCL4","CCL5", 
                                                             "PDCD1","CTLA4","TIGIT","FOXP3", "CD8A"))

pdf("SCS_SNS_CD4_heatmap_selected.pdf", width = 10, height = 10)
ggplot(avgexp.heatmap, aes(x = gene, y = cell.type, fill = avg.expression)) +
  geom_tile(color = "white") +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradientn(colors = c( "#f7fbff", "#deebf7", "#4292c6", "#2171b5", "#08519c", "#08306b")) 
dev.off()


# Fig 2E -----------------------------------------------------------------------

ggplot1 <- DimPlot(CD8.combined, label = FALSE, group.by = "annotation", pt.size=1,
                   cols = c('TEM-a' = '#bf8fc9',
                            
                            'TEM-b' = '#9ecae1',
                            'TEM-c' = '#807dba',
                            'Naive' = '#5a1209',
                            'TEM-d' = '#ed8f87',
                            'TCM-b' = '#8c564b',
                            
                            'TEM-e' = '#f5af2d',
                            'MT-high' = 'lightgrey',
                            'TEM-f' = '#bcbd22',
                            'CD4' = '#a50f15',
                            'TCM-a' = '#225ea8',
                            'TEM-g' = '#f3bd91'))
ggsave("SCS_SNS_CD8_umap.jpg", plot = ggplot1, width = 10, height = 10, dpi = 300)

DefaultAssay(CD8.combined) <- "RNA"
genes.heatmap <- c("S1PR1", "CCR7", "SELL", "FAS",
                   "CXCR3", "CCR5", "ITGAL",
                   "HLA-DRA", "HLA-DRB1", "HLA-DQB1",
                   "FOS","JUN", "CD40LG",
                   "CCL4","CCL5",
                   "PDCD1","CTLA4","TIGIT",
                   "MALAT1", "MT-CO1", "CD4")
CD8.combined <- ScaleData(CD8.combined,features = genes.heatmap)

avgexp <- as.data.frame(AverageExpression(CD8.combined,features = genes.heatmap,assays = "RNA",slot = "scale.data")$RNA)
avgexp.heatmap <- avgexp[rownames(avgexp) %in% genes.heatmap,] %>% rownames_to_column(var = "gene") %>% gather(-gene,key = "cell.type",value = "avg.expression")

avgexp.heatmap$cell.type = factor(avgexp.heatmap$cell.type, levels = rev(c("Naive", "TCM-a", "TCM-b", "TEM-a", "TEM-b", "TEM-c", "TEM-d", "TEM-e", "TEM-f", "TEM-g", "MT-high", "CD4")))
avgexp.heatmap$gene = factor(avgexp.heatmap$gene, levels = c("S1PR1", "CCR7", "SELL", "FAS",
                                                             "CXCR3", "CCR5", "ITGAL",
                                                             "HLA-DRA", "HLA-DRB1", "HLA-DQB1",
                                                             "FOS","JUN", "CD40LG",
                                                             "CCL4","CCL5",
                                                             "PDCD1","CTLA4","TIGIT",
                                                             "MALAT1", "MT-CO1", "CD4"))

pdf("SCS_SNS_CD8_heatmap_selected.pdf", width = 10, height = 10)
ggplot(avgexp.heatmap, aes(x = gene, y = cell.type, fill = avg.expression)) +
  geom_tile(color = "white") +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradientn(colors = c("#fcfbfd", "#efedf5",  "#807dba", "#6a51a3", "#54278f", "#3f007d"))
dev.off()

# Fig 2F -----------------------------------------------------------------------

ggplot1 <- DimPlot(hepatocytes.combined, label = FALSE, group.by = "annotation", pt.size=1,
                   cols = c('Central' = '#119891',
                            'Portal' = '#005679',
                            'Midzone' = '#8770ac',
                            'Inflammatory' = '#1e8adb',
                            'BICC1-hi' = '#65a644',
                            'Proliferation-2' = '#c87726',
                            'Proliferation-1' = '#9f6763',
                            
                            'PTPRC-hi' = '#969696',
                            'PTPRB-hi' = '#525252',
                            'Un-define' = 'lightgrey'))
ggsave("SNS_hepatocytes_umap.jpg", plot = ggplot1, width = 10, height = 10, dpi = 300)

DefaultAssay(hepatocytes.combined) <- "RNA"
genes.heatmap <- c("SDS", "HAL", "PCK1", 
                   "CYP1A2","CYP2E1",
                   "GLUL","LGR5", 
                   "IL32", "ICAM1", "STAT1",
                   "DIAPH3","BRIP1",                   
                   "MKI67",
                   "BICC1","DCDC2", "FGF13",
                   "PTPRB", "PTPRC")
hepatocytes.combined <- ScaleData(hepatocytes.combined,features = genes.heatmap)

avgexp <- as.data.frame(AverageExpression(hepatocytes.combined,features = genes.heatmap,assays = "RNA",slot = "scale.data")$RNA)
head(avgexp)

avgexp.heatmap <- avgexp[rownames(avgexp) %in% genes.heatmap,] %>% rownames_to_column(var = "gene") %>% gather(-gene,key = "cell.type",value = "avg.expression")
head(avgexp.heatmap)

avgexp.heatmap$cell.type = factor(avgexp.heatmap$cell.type, levels = rev(c("Portal", "Midzone", "Central", "Inflammatory", "Proliferation-1", "Proliferation-2", "BICC1-hi", "PTPRB-hi", "PTPRC-hi", "Un-define")))
avgexp.heatmap$gene = factor(avgexp.heatmap$gene, levels = c( "SDS", "HAL", "PCK1", 
                                                              "CYP1A2","CYP2E1",
                                                              "GLUL","LGR5", 
                                                              "IL32", "ICAM1", "STAT1",
                                                              "DIAPH3","BRIP1",                   
                                                              "MKI67",
                                                              "BICC1","DCDC2", "FGF13",
                                                              "PTPRB", "PTPRC"))

pdf("SNS_hepatocytes_heatmap_selected.pdf", width = 10, height = 10)
ggplot(avgexp.heatmap, aes(x = gene, y = cell.type, fill = avg.expression)) +
  geom_tile(color = "white") +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradientn(colors = c("#f7fcf5", "#e5f5e0",  "#41ab5d", "#238b45", "#006d2c", "#00441b"))
dev.off()