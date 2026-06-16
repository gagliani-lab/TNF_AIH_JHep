## Extended Figure 3 ##

library(Seurat)
library(dplyr)
library(ggplot2)

hepatocytes.combined <- readRDS(file='./Integration/RDS/14_SNS_hepatocytes_figures.rds')
CD4.combined <- readRDS(file='./Integration/RDS/12_SCS_SNS_CD4_figures.rds')
CD8.combined <- readRDS(file='./Integration/RDS/13_SCS_SNS_CD8_figures.rds')

# ED Fig 3A --------------------------------------------------------------------

DefaultAssay (CD4.combined) <- "RNA"

cluster_color <- c('lightgrey', #Naive
                   '#d99f3e', #TCM
                   '#c18f6e', #TEM-a
                   '#ae98b6', #TEM-b
                   '#003c96', #TRM1
                   '#41ab5d', #TFH
                   '#426e66', #TR1
                   '#807dba', #TREG
                   '#dadaeb' #TREG-naive
)

ggplot2 <- DoHeatmap(CD4.combined, 
                     features = c("CCR7", "ACTN1", "SELL", "TCF7", "LEF1", "LYPD3", "EEF1B2", #Naive
                                  "IL7R", "ANXA1", "LMNA", "MAL", "KLF2", "S1PR1", "CD55", #TCM
                                  "HLA-DRB5", "HLA-DQB1", "HLA-DRA", "HLA-DQA1", "HOPX", "HLA-DRB1", "HLA-DPA1", #TEM
                                  "KLRK1", "AOAH", "PECAM1", "NKG7", "ZEB2", "CCL5", "ARHGAP26", #TEM-NKG7-high
                                  "CCL5", "CD69", "GZMA", "IFNG", "KLRB1", "TNF", "KLF6",#TRM1
                                  "TOX2", "RNF19A", "ST8SIA1", "NR3C1", "TOX", "PDCD1", "TNFSF8",#TFH
                                  "IL10", "GZMK", "ENC1", "HAVCR2", "CD38", "DTHD1", "GZMA",#TR1
                                  "FOXP3", "IKZF2", "TNFRSF9", "CCR8", "LAYN", "AC017002.3", "TBC1D4",#TREG
                                  "FOXP3", "RTKN2", "IKZF2", "CTLA4", "TIGIT", "STAM", "PLCL1" #TREG-naive                      
                     ),
                     size = 4, 
                     slot = "scale.data", 
                     group.colors = cluster_color, 
                     draw.lines = TRUE,
                     lines.width = 100)
ggsave("SCS_SNS_CD4_heatmap_DEG.jpg", plot = ggplot2, width = 10, height = 10, dpi = 300)


# ED Fig 3B --------------------------------------------------------------------

proportions_2 <- as.data.frame(prop.table(table(immune.combined$annotation, immune.combined$orig.ident),margin = 1))
colnames(proportions_2) <- c("annotation","sample","proportion")
proportions_2$annotation <- factor(proportions_2$annotation, 
                                   levels = c("Naive", "TCM", "TEM-a", "TEM-b", "TRM1", "TFH", "TR1", "TREG", "TREG-naive"))

pdf("SCS_SNS_CD4_proportion_annotation.pdf", width = 10, height = 10)
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


proportions_1 <- as.data.frame(prop.table(table(immune.combined$annotation, immune.combined$orig.ident),margin = 2))
colnames(proportions_1) <- c("annotation","sample","proportion")
proportions_1$annotation <- factor(proportions_1$annotation, 
                                   levels = c("Naive", "TCM", "TEM-a", "TEM-b", "TRM1", "TFH", "TR1", "TREG", "TREG-naive"))
pdf("SCS_SNS_CD4_proportion_sample.pdf", width = 10, height = 10)
ggplot(proportions_1, aes(x = sample, y = proportion, fill = annotation)) +
        geom_col() +
        theme_classic() +
        scale_fill_manual(values = c('lightgrey', #Naive
                                     '#d99f3e', #TCM
                                     '#c18f6e', #TEM-a
                                     '#ae98b6', #TEM-b
                                     '#003c96', #TRM1
                                     '#41ab5d', #TFH
                                     '#426e66', #TR1
                                     '#807dba', #TREG
                                     '#dadaeb' #TREG-naive 
        ))
dev.off()

# ED Fig 3C --------------------------------------------------------------------
DefaultAssay (CD8.combined) <- "RNA"

cluster_color <- c('#5a1209', #Naive
                   '#225ea8', #TCM-a
                   '#8c564b', #TCM-b
                   '#bf8fc9', #TEM-a
                   '#9ecae1', #TEM-b
                   '#807dba', #TEM-c
                   '#ed8f87', #TEM-d
                   '#f5af2d', #TEM-e
                   '#bcbd22', #TEM-f
                   '#f3bd91', #TEM-g
                   'lightgrey', #MT-high
                   '#a50f15' #CD4
)
ggplot2 <- DoHeatmap(CD8.combined, 
                     features = c("SELL", "CCR7", "LEF1", "ACTN1", "MAL", "LYPD3", "CD248", #Naive
                                  "GNLY", "FGFBP2", "GZMB", "FCGR3A", "CX3CR1", "LGALS1", "EFHD2", #TCM-a
                                  "IL7R", "LTB", "GPR183", "TCF7", "LEF1", "SELL", "ANXA1", #TCM-b
                                  "CMC1", "GZMK", "CD160", "IFNGR1", "KLRG1", "GZMM", "AUTS2", #TEM-a
                                  "FCER1G", "TYROBP", "NCAM1", "B3GNT7", "TMIGD2", "KLRF1", "SH2D1B", #TEM-b
                                  "TXNIP", "CCL5", "GIMAP7", "PTPRCAP", "MT-CO1", "GIMAP4", "HCST", #TEM-c
                                  "THEMIS", "AOAH", "RUNX1", "LINC01934", "SLFN12L", "CHST11", #TEM-d 
                                  "CXCR6", "CTLA4", "SIRPG", "HAVCR2", "LAYN", "TNFSF4", "RGS1", #TEM-e
                                  "CCL4L2", "CCL4", "IFNG", "TNF", "TNFSF9", "CCL3", "FOS", #TEM-f 
                                  "DGKH", "SNX9", "GFOD1", "FAM3C", "LYST", "TNFRSF9", "PTPN11", #TEM-g 
                                  "PTPRC", "RNF213", "TNFAIP3", "MALAT1", "MT-CO1", "MACF1", "SLC38A1", #MT-high
                                  "FOXP3", "TNFRSF18", "CCR8", "IL2RA", "IL1R2", "CD4", "TNFRSF4" #FOXP3+                       
                     ),
                     size = 4, 
                     slot = "scale.data", 
                     group.colors = cluster_color, 
                     draw.lines = TRUE,
                     lines.width = 100)
ggsave("SCS_SNS_CD8_heatmap_DEG.jpg", plot = ggplot2, width = 10, height = 10, dpi = 300)

# ED Fig 3D --------------------------------------------------------------------

proportions_2 <- as.data.frame(prop.table(table(CD8.combined$annotation, CD8.combined$orig.ident),margin = 1))
colnames(proportions_2) <- c("annotation","sample","proportion")
# Reordering the 'annotation' column based on a specific order
proportions_2$annotation <- factor(proportions_2$annotation, 
                                   levels = c("Naive", "TCM-a", "TCM-b", "TEM-a", "TEM-b", "TEM-c", "TEM-d", "TEM-e", "TEM-f", "TEM-g", "MT-high", "CD4"))
pdf("SCS_SNS_CD8_proportion_annotation.pdf", width = 10, height = 10)
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

proportions_1 <- as.data.frame(prop.table(table(CD8.combined$annotation, CD8.combined$orig.ident),margin = 2))
head(proportions_1)

colnames(proportions_1) <- c("annotation","sample","proportion")
# Reordering the 'annotation' column based on a specific order
proportions_1$annotation <- factor(proportions_1$annotation, 
                                   levels = c("Naive", "TCM-a", "TCM-b", "TEM-a", "TEM-b", "TEM-c", "TEM-d", "TEM-e", "TEM-f", "TEM-g", "MT-high", "CD4"))

pdf("SCS_SNS_CD8_proportion_sample.pdf", width = 10, height = 10)
ggplot(proportions_1, aes(x = sample, y = proportion, fill = annotation)) +
        geom_col() +
        theme_classic() +
        scale_fill_manual(values = c('#5a1209', #Naive
                                     '#225ea8', #TCM-a
                                     '#8c564b', #TCM-b
                                     '#bf8fc9', #TEM-a
                                     '#9ecae1', #TEM-b
                                     '#807dba', #TEM-c
                                     '#ed8f87', #TEM-d
                                     '#f5af2d', #TEM-e
                                     '#bcbd22', #TEM-f
                                     '#f3bd91', #TEM-g
                                     'lightgrey', #MT-high
                                     '#a50f15' #CD4
        ))
dev.off()


# ED Fig 3E --------------------------------------------------------------------

DefaultAssay (hepatocytes.combined) <- "RNA"

cluster_color <- c('#005679', #Portal
                   '#8770ac', #Midzone
                   '#119891', #Central
                   '#1e8adb', #Inflammation
                   '#9f6763', #Proliferation-1
                   '#c87726', #Proliferation-2
                   '#65a644', #BICC1-hi
                   '#525252', #PTPRB-hi
                   '#969696', #PTPRC-hi
                   'lightgrey' #Un-define                 
)

ggplot2 <- DoHeatmap(hepatocytes.combined, 
                     features = c("SDS", "HAL", "CYP2A7", "ASS1", "LINC00598", "UPP2", "CCNT2-AS1", #Portal
                                  "CYP2C8", "CYP2C9", "TF", "EBNA1BP2", "ALB", "PLG", "ADH1B", #Midzone
                                  "RELN", "LGR5", "GLUL", "SLCO1B3", "TENM2", "LINC01344", "RHBG", #Central
                                  "CXCL10", "AC083837.1", "GBP1", "BIRC3", "GBP4", "OR2I1P", "IL32", "ICAM1", #Inflammatory
                                  "MIR924HG", "DIAPH3", "ASPM", "BRIP1", "CENPP", "APOLD1", "EZH2", #Proliferative-1
                                  "CENPF", "HIST1H2BD", "HIST1H2BC", "MKI67", "TOP2A", "AL353759.1", "CENPE", #Proliferative-2
                                  "BICC1", "CTNND2", "NRSN1", "CDH6", "LINC02331", "WNK2", "FGFR2", "DCDC2", #BICC1-hi
                                  "LDB2", "FBXL7", "ST6GALNAC3", "PTPRB", "BMPER", "DNASE1L3", "CRHBP", #PTPRB-hi
                                  "PTPRC", "ARHGAP15", "CELF2", "TOX", "AC079793.1", "PRKCH", "FYN", #PTPRC-hi
                                  "HOOK3", "AMBP", "ARHGEF7" , "ANKS1A", "DDX46", "HP1BP3", "NUFIP2" #Un-define
                     ),
                     size = 4, 
                     slot = "scale.data", 
                     group.colors = cluster_color, 
                     draw.lines = TRUE,
                     lines.width = 100)
ggsave("SNS_hepatocytes_heatmap_DEG.jpg", plot = ggplot2, width = 10, height = 10, dpi = 300)

# ED Fig 3F --------------------------------------------------------------------

proportions_2 <- as.data.frame(prop.table(table(hepatocytes.combined$annotation, hepatocytes.combined$orig.ident),margin = 1))
colnames(proportions_2) <- c("annotation","sample","proportion")
proportions_2$annotation <- factor(proportions_2$annotation, 
                                   levels = c("Portal", "Midzone", "Central", "Inflammatory", "Proliferation-1", "Proliferation-2", "BICC1-hi", "PTPRB-hi", "PTPRC-hi", "Un-define"))

pdf("SNS_hepatocytes_proportion_annotation.pdf", width = 10, height = 10)
ggplot(proportions_2, aes(x=annotation,y=proportion,fill=sample)) +
        geom_col() +
        theme_classic() + scale_fill_manual(values = c(
                'NUC0008' = '#08306b',
                'NUC010A' = '#08519c', 
                'NUC015A' = '#2171b5',
                'NUC0240' = '#4292c6',
                'NUC0691' = '#6baed6',
                'NUC0786' = '#9ecae1'))
dev.off()

proportions_1 <- as.data.frame(prop.table(table(hepatocytes.combined$annotation, hepatocytes.combined$orig.ident),margin = 2))
colnames(proportions_1) <- c("annotation","sample","proportion")
proportions_1$annotation <- factor(proportions_1$annotation, 
                                   levels = c("Portal", "Midzone", "Central", "Inflammatory", "Proliferation-1", "Proliferation-2", "BICC1-hi", "PTPRB-hi", "PTPRC-hi", "Un-define"))

pdf("SNS_hepatocztes_proportion_sample.pdf", width = 10, height = 10)
ggplot(proportions_1, aes(x = sample, y = proportion, fill = annotation)) +
        geom_col() +
        theme_classic() +
        scale_fill_manual(values = c('#005679', #Portal
                                     '#8770ac', #Midzone
                                     '#119891', #Central
                                     '#1e8adb', #Inflammation
                                     '#9f6763', #Proliferation-1
                                     '#c87726', #Proliferation-2
                                     '#65a644', #BICC1-hi
                                     '#525252', #PTPRB-hi
                                     '#969696', #PTPRC-hi
                                     'lightgrey' #Un-define
        ))
dev.off()