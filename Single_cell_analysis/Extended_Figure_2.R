## Extended Figure 2 ##

library(Seurat)
library(dplyr)
library(ggplot2)

AIH.combined <- readRDS(file='./Integration/RDS/3_SCS_SNS_atlas.rds')
SCS.combined <- readRDS(file='./Integration/RDS/5_SCS_atlas.rds')
SNS.combined <- readRDS(file='./Integration/RDS/6_SNS_atlas.rds')

# ED Fig 2A --------------------------------------------------------------------


AIH.combined$nCount.capped <- ifelse(AIH.combined$nCount_RNA>40000,40000,AIH.combined$nCount_RNA)
AIH.combined$nFeature.capped <- ifelse(AIH.combined$nFeature_RNA>10000,10000,AIH.combined$nFeature_RNA)

pdf("atlas_nFeature.pdf", width = 10, height = 10)
VlnPlot(AIH.combined, features = c("nFeature.capped"), 
        split.by = "protocol", 
        cols= c('#238b45','#cc338b'), 
        pt.size=0.01,
        split.plot = TRUE) 
dev.off()

pdf("atlas_nCount.pdf", width = 10, height = 10)
VlnPlot(AIH.combined, features = c("nCount.capped"), 
        split.by = "protocol",
        cols= c('#238b45','#cc338b'), 
        pt.size=0.01,
        split.plot = TRUE)
dev.off()

ggplot_nFeature <- VlnPlot(AIH.combined, features = c("nFeature.capped"), 
                           split.by = "protocol", 
                           cols= c('#238b45','#cc338b'), 
                           pt.size=0.01,
                           split.plot = TRUE) 
ggsave("atlas_nFeature.jpg", plot = ggplot_nFeature, width = 10, height = 10, dpi = 300)

ggplot_nCount <- VlnPlot(AIH.combined, features = c("nCount.capped"), 
                         split.by = "protocol", 
                         cols= c('#238b45','#cc338b'), 
                         pt.size=0.01,
                         split.plot = TRUE) 
ggsave("atlas_nCount.jpg", plot = ggplot_nCount, width = 10, height = 10, dpi = 300)

# ED Fig 2C --------------------------------------------------------------------

ggplot_1 <- DimPlot(SCS.combined, reduction = "umap", label = FALSE,
                    cols = c('CD8' = '#6baed6',
                             'CD4' = '#08519c',
                             'Myeloid' = '#df65b0',
                             'Innate_like_T' = '#1d91c0',
                             'Innate_2' = '#EF6538',
                             'Innate_1' = '#EF6538',
                             'B' = '#fd8d3c',
                             'Plasma_cell' = '#cb181d',
                             'Proliferating' = '#c6dbef',
                             'MT-hi' = '#737373'))
ggsave("SCS_atlas_umap.jpg", plot = ggplot_1, width = 10, height = 10, dpi = 300)

projection_SCS <- read.csv("./Integration/Table/SCS_atlas_projection_annotation.csv")
rownames(projection_SCS) <- projection_SCS$X
projection_SCS$X <- NULL
test_SNS_SCS <- AIH.combined
test_SNS_SCS <- AddMetaData(object = test_SNS_SCS, metadata = projection_SCS, col.name= colnames(projection_SCS))
test_SNS_SCS$SCS_annotation[is.na(test_SNS_SCS$SCS_annotation)] <- "SNS"
Idents(test_SNS_SCS) <- "SCS_annotation"

ggplot_4 <- DimPlot(test_SNS_SCS, label=FALSE, group.by = "SCS_annotation", split.by = "protocol",
                    cols = c('CD8' = '#6baed6',
                             'CD4' = '#08519c',
                             'Myeloid' = '#df65b0',
                             'Innate_like_T' = '#1d91c0',
                             'Innate_2' = '#EF6538',
                             'Innate_1' = '#EF6538',
                             'B' = '#fd8d3c',
                             'Plasma_cell' = '#cb181d',
                             'Proliferating' = '#c6dbef',
                             'MT-hi' = '#737373'))

ggsave("Atlas_SCS_projection_umap.jpg", plot = ggplot_4, width = 17, height = 10, dpi = 300)

# ED Fig 2D --------------------------------------------------------------------

ggplot_1 <- DimPlot(SNS.combined, reduction = "umap", label = FALSE,
                    cols = c('CD8' = '#6baed6',
                             'Hepatocytes' = '#006837',
                             'Proliferating' = '#c6dbef',
                             'CD4' = '#08519c',
                             'Myeloid' = '#df65b0',
                             'Innate' = '#EF6538',
                             'Endothelial' = '#addd8e',
                             'Cholangiocytes' = '#41ab5d',
                             'Stellate/Fibroblasts' = '#807dba', 
                             'BICC1-hi_hepato' = '#65a644', 
                             'B' = '#fd8d3c',
                             'Plasma' = '#cb181d', 
                             'Undefine_1' = '#737373'))
ggsave("SNS_atlas_umap.jpg", plot = ggplot_1, width = 10, height = 10, dpi = 300)

projection_SNS <- read.csv("C:/Yang_AIH_Final_Version/Integration/Table/SNS_atlas_projection_annotation.csv")
rownames(projection_SNS) <- projection_SNS$X
projection_SNS$X <- NULL
test_SNS_SCS <- AIH.combined
test_SNS_SCS <- AddMetaData(object = test_SNS_SCS, metadata = projection_SNS, col.name= colnames(projection_SNS))
test_SNS_SCS$SNS_annotation[is.na(test_SNS_SCS$SNS_annotation)] <- "SCS"
Idents(test_SNS_SCS) <- "SNS_annotation"

ggplot_3 <- DimPlot(test_SNS_SCS, label=FALSE, group.by = "SNS_annotation", split.by = "protocol",
                    cols = c('CD8' = '#6baed6',
                             'Hepatocytes' = '#006837',
                             'Proliferating' = '#c6dbef',
                             'CD4' = '#08519c',
                             'Myeloid' = '#df65b0',
                             'Innate' = '#EF6538',
                             'Endothelial' = '#addd8e',
                             'Cholangiocytes' = '#41ab5d',
                             'Stellate/Fibroblasts' = '#807dba', 
                             'BICC1-hi_hepato' = '#65a644', 
                             'B' = '#fd8d3c',
                             'Plasma' = '#cb181d', 
                             'Undefine_1' = '#737373'))
ggsave("Atlas_SNS_projection_umap.jpg", plot = ggplot_3, width = 17, height = 10, dpi = 300)


# ED Fig 2E --------------------------------------------------------------------

cluster_color <- c('#6baed6', #CD8
                   '#08519c', #CD4
                   '#df65b0', #Myeloid
                   '#EF6538', #Innate_1
                   '#EF6538', #Innate_2
                   '#1d91c0', #Innate_like_T
                   '#fd8d3c', #B
                   '#cb181d', #Plasma_cell
                   '#006837', #Hepatocytes
                   '#addd8e', #Endothelial
                   '#41ab5d', #Cholangiocytes
                   '#807dba', #Stellate/Fibroblasts
                   '#c6dbef', #Proliferating
                   '#737373', #MT-high
                   '#737373', #Undefine_1
                   '#737373'  #Undefine_2
)
ggplot3 <- DoHeatmap(AIH.combined, 
                     features = c("CD8A", "CCL5", "TNFSF9", "CD8B", "RGS1", #CD8
                                  "IL7R","LTB","CD52","ANXA1","KLF2", #CD4
                                  "FMNL2","SLC8A1","KCNMA1","CD163","IFI30", #Myeloid
                                  "XCL1","CCL3","AREG","FCER1G","TYROBP", #Innate_1
                                  "GNLY","GZMB","FGFBP2","PRF1","FCGR3A", #Innate_2
                                  "IL7R","SLC4A10","LTB","KLRB1","CD69", #Innate_like_T
                                  "MS4A1", "CD79A", "HLA-DRA", "CD74", "BANK1", #B
                                  "IFNG-AS1", "IGHG1", "IGKC", "TXNDC5", "PDK1", #Plasma_cell
                                  "CPS1", "ACSM2A", "SLCO1B3", "TF", "CYP2C9", #Hepatocytes
                                  "LDB2", "NRG3", "LIFR", "PTPRB", "STAB2", #Endothelial
                                  "BICC1", "PKHD1", "CFTR", "ANXA4", "FGFR2", #Cholangiocytes
                                  "RYR2", "PRKG1", "LAMA2", "CACNA1C", "CCBE1", #Stellate/Fibroblasts
                                  "STMN1", "RRM2", "HMGB2", "MKI67", "ASPM", #Proliferating
                                  "TNFAIP3", "PTPRC", "MALAT1", "MT-CO1", "RNF19A", #MT-high
                                  "SLFN12L", "PRKCH", "TOX", "AOAH", "MBNL1", #Undefine_1 
                                  "PLD4", "LINC01478", "JCHAIN", "BCL11A", "LILRA4"
                     ),
                     size = 4, 
                     slot = "scale.data", 
                     group.colors = cluster_color, 
                     draw.lines = TRUE,
                     lines.width = 100)
ggsave("SCS_SNS_atlas_DEG_heatmap_deg.jpg", plot = ggplot3, width = 10, height = 10, dpi = 300)

# ED Fig 2F --------------------------------------------------------------------

proportions_2 <- as.data.frame(prop.table(table(AIH.combined$annotation, AIH.combined$orig.ident),margin = 1))
colnames(proportions_2) <- c("annotation","sample","proportion")

proportions_2$annotation <- factor(proportions_2$annotation, 
                                   levels = c('CD8', 'CD4', 'Myeloid', 'Innate_1', 'Innate_2', 
                                              'Innate_like_T', 'B', 'Plasma_cell', 'Hepatocytes', 
                                              'Endothelial', 'Cholangiocytes', 'Stellate/Fibroblasts', 
                                              'Proliferating', 'MT-high', 'Undefine_1', 'Undefine_2'))

pdf("SCS_SNS_proportion_annotation.pdf", width = 10, height = 10)
ggplot(proportions_2, aes(x=annotation,y=proportion,fill=sample)) +
        geom_col() +
        theme_classic() + scale_fill_manual(values = c('A05' = '#efedf5',
                                                       'A07' = '#dadaeb',
                                                       'A08' = '#bcbddc',
                                                       'A09' = '#9e9ac8',
                                                       'A10' = '#807dba', 
                                                       'A11' = '#6a51a3',
                                                       'A12' = '#54278f',
                                                       'A13' = '#3f007d',
                                                       'A14' = '#810f7c',
                                                       'A15' = '#88419d',
                                                       'NUC0008' = '#08306b',
                                                       'NUC010A' = '#08519c', 
                                                       'NUC015A' = '#2171b5',
                                                       'NUC0240' = '#4292c6',
                                                       'NUC0691' = '#6baed6',
                                                       'NUC0786' = '#9ecae1'))
dev.off()

proportions_1 <- as.data.frame(prop.table(table(AIH.combined$annotation, AIH.combined$orig.ident),margin = 2))
colnames(proportions_1) <- c("annotation","sample","proportion")
# Reordering the 'annotation' column based on a specific order
proportions_1$annotation <- factor(proportions_1$annotation, 
                                   levels = c('CD8', 'CD4', 'Myeloid', 'Innate_1', 'Innate_2', 
                                              'Innate_like_T', 'B', 'Plasma_cell', 'Hepatocytes', 
                                              'Endothelial', 'Cholangiocytes', 'Stellate/Fibroblasts', 
                                              'Proliferating', 'MT-high', 'Undefine_1', 'Undefine_2'))
pdf("SCS_SNS_proportion_sample.pdf", width = 10, height = 10)
ggplot(proportions_1, aes(x = sample, y = proportion, fill = annotation)) +
        geom_col() +
        theme_classic() +
        scale_fill_manual(values = c('CD8' = '#6baed6',
                                     'CD4' = '#08519c',
                                     'Myeloid' = '#df65b0',
                                     'Innate_1' = '#EF6538',
                                     'Innate_2' = '#EF6538', 
                                     'Innate_like_T' = '#1d91c0',
                                     'B' = '#fd8d3c',
                                     'Plasma_cell' = '#cb181d',
                                     'Hepatocytes' = '#006837',
                                     'Endothelial' = '#addd8e',
                                     'Cholangiocytes' = '#41ab5d',
                                     'Stellate/Fibroblasts' = '#807dba', 
                                     'Proliferating' = '#c6dbef',
                                     'MT-high' = '#737373',
                                     'Undefine_1' = '#737373',
                                     'Undefine_2' = '#737373'))
dev.off()