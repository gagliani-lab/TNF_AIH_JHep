## Integration + Annotation AIH atlas ##

library(Seurat)
library(dplyr)

SCS_05 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH05.rds")
SCS_05$protocol <- "SCS"
SCS_07 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH07.rds")
SCS_07$protocol <- "SCS"
SCS_08 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH08.rds")
SCS_08$protocol <- "SCS"
SCS_09 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH09.rds")
SCS_09$protocol <- "SCS"
SCS_10 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH10.rds")
SCS_10$protocol <- "SCS"
SCS_11 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH11.rds")
SCS_11$protocol <- "SCS"
SCS_12 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH12.rds")
SCS_12$protocol <- "SCS"
SCS_13 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH13.rds")
SCS_13$protocol <- "SCS"
SCS_14 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH14.rds")
SCS_14$protocol <- "SCS"
SCS_15 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_AIH15.rds")
SCS_15$protocol <- "SCS"

SNS_0008 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_NUC0008.rds")
SNS_0008$protocol <- "SNS"
SNS_010A <- readRDS("./Single_sample_preprocess/RDS/prefiltered_NUC010A.rds")
SNS_010A$protocol <- "SNS"
SNS_015A <- readRDS("./Single_sample_preprocess/RDS/prefiltered_NUC015A.rds")
SNS_015A$protocol <- "SNS"
SNS_0240 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_NUC0240.rds")
SNS_0240$protocol <- "SNS"
SNS_0691 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_NUC0691.rds")
SNS_0691$protocol <- "SNS"
SNS_0786 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_NUC0786.rds")
SNS_0786$protocol <- "SNS"

merge <- merge(SCS_05, y = c(SCS_07, 
                             SCS_08, 
                             SCS_09, 
                             SCS_10, 
                             SCS_11, 
                             SCS_12, 
                             SCS_13, 
                             SCS_14, 
                             SCS_15, 
                             SNS_0008, 
                             SNS_010A, 
                             SNS_015A, 
                             SNS_0240, 
                             SNS_0691, 
                             SNS_0786), 
               add.cell.ids = c("AIH05", 
                                "AIH07", 
                                "AIH08", 
                                "AIH09", 
                                "AIH10", 
                                "AIH11", 
                                "AIH12", 
                                "AIH13", 
                                "AIH14", 
                                "AIH15",
                                "NUC0008",
                                "NUC010A",
                                "NUC015A",
                                "NUC0240",
                                "NUC0691",
                                "NUC0786"))

AIH05.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH05_singlet_barcodes.txt",header = FALSE)$V1
AIH07.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH07_singlet_barcodes.txt",header = FALSE)$V1
AIH08.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH08_singlet_barcodes.txt",header = FALSE)$V1
AIH09.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH09_singlet_barcodes.txt",header = FALSE)$V1
AIH10.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH10_singlet_barcodes.txt",header = FALSE)$V1
AIH11.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH11_singlet_barcodes.txt",header = FALSE)$V1
AIH12.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH12_singlet_barcodes.txt",header = FALSE)$V1
AIH13.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH13_singlet_barcodes.txt",header = FALSE)$V1
AIH14.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH14_singlet_barcodes.txt",header = FALSE)$V1
AIH15.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/AIH15_singlet_barcodes.txt",header = FALSE)$V1
SNS0008.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0008_singlet_barcodes.txt",header = FALSE)$V1
SNS010A.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC010A_singlet_barcodes.txt",header = FALSE)$V1
SNS015A.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC015A_singlet_barcodes.txt",header = FALSE)$V1
SNS0240.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0240_singlet_barcodes.txt",header = FALSE)$V1
SNS0691.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0691_singlet_barcodes.txt",header = FALSE)$V1
SNS0786.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0786_singlet_barcodes.txt",header = FALSE)$V1

all.SCS.singlet.barcodes <- c(AIH05.singlet.barcodes, 
                              AIH07.singlet.barcodes,
                              AIH08.singlet.barcodes,
                              AIH09.singlet.barcodes,
                              AIH10.singlet.barcodes,
                              AIH11.singlet.barcodes,
                              AIH12.singlet.barcodes,
                              AIH13.singlet.barcodes,
                              AIH14.singlet.barcodes,
                              AIH15.singlet.barcodes,
                              SNS0008.singlet.barcodes,
                              SNS010A.singlet.barcodes,
                              SNS015A.singlet.barcodes,
                              SNS0240.singlet.barcodes,
                              SNS0691.singlet.barcodes,
                              SNS0786.singlet.barcodes)

merge$DF.status <- ifelse(rownames(merge@meta.data) %in% all.SCS.singlet.barcodes,"Singlet","Doublet")

merge_singlets <- subset(merge, subset = DF.status == "Singlet")

sample.list <- SplitObject(merge, split.by = "orig.ident")
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = sample.list)
AIH.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
AIH.combined <- IntegrateData(anchorset = AIH.anchors)

AIH.combined <- readRDS("C://Yang_AIH_Final_Version//Integration//RDS//2_SCS_SNS_FindAnchor_step1.rds")

DefaultAssay(AIH.combined) <- "integrated"
AIH.combined <- ScaleData(AIH.combined, verbose = FALSE)
AIH.combined <- RunPCA(AIH.combined, npcs = 30, verbose = FALSE)
AIH.combined <- RunUMAP(AIH.combined, reduction = "pca", dims = 1:30)
AIH.combined <- FindNeighbors(AIH.combined, reduction = "pca", dims = 1:30)
AIH.combined <- FindClusters(AIH.combined, resolution = 0.8)

DefaultAssay(AIH.combined) <- "RNA"

AIH.combined <- NormalizeData(AIH.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
AIH.combined <- ScaleData(AIH.combined)

DefaultAssay(AIH.combined) <- "CITE"

AIH.combined <- NormalizeData(AIH.combined, normalization.method = "CLR", margin = 2, assay = "CITE")
AIH.combined <- ScaleData(AIH.combined)


new.cluster.ids <- c(
  'CD8', #0
  
  'Hepatocytes', #1
  'CD4', #2
  'CD8', #3
  'Innate_1', #4
  'Endothelial', #5
  
  'Myeloid', #6
  'Hepatocytes', #7
  'Hepatocytes', #8
  'Myeloid', #9
  'Undefine_1', #10
  
  'Cholangiocytes', #11
  'Innate_2', #12
  'CD8', #13
  'Hepatocytes', #14
  'Innate_like_T', #15
  'Endothelial', #16
  'Hepatocytes', #17
  'Stellate/Fibroblasts', #18
  'B', #19
  'CD4', #20
  'Proliferating', #21
  'Plasma_cell', #22
  'MT-high', #23
  'Myeloid',#24
  'Undefine_2' #25
)

names(new.cluster.ids) <- levels(AIH.combined)
AIH.combined <- RenameIdents (AIH.combined, new.cluster.ids)
AIH.combined$annotation <- AIH.combined@active.ident

write.table(data.frame(barcode = colnames(AIH.combined),celltype = AIH.combined$annotation),"./Integration/barcodes_SCS_SNS_atlas.txt",sep="\t",row.names = FALSE,quote = FALSE)

AIH.combined.markers <- FindAllMarkers(AIH.combined, assay = "RNA",logfc.threshold = 0.5,min.pct = 0.25,only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(AIH.combined.markers, file = "./Integration/DEG_table/SCS_SNS_atlas_DEG_cluster.txt", sep = "\t",quote=FALSE, row.names=FALSE)

saveRDS(AIH.combined, file='./Integration/RDS/3_SCS_SNS_atlas.rds')

