## Integration + Annotation SCS atlas ##

library(dplyr)
library(Seurat)
library(ggplot2)

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

merge <- merge(SNS_0008, y = c(SNS_010A, 
                               SNS_015A, 
                               SNS_0240, 
                               SNS_0691, 
                               SNS_0786), 
               add.cell.ids = c("NUC0008",
                                "NUC010A",
                                "NUC015A",
                                "NUC0240",
                                "NUC0691",
                                "NUC0786"))

SNS0008.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0008_singlet_barcodes.txt",header = FALSE)$V1
SNS010A.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC010A_singlet_barcodes.txt",header = FALSE)$V1
SNS015A.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC015A_singlet_barcodes.txt",header = FALSE)$V1
SNS0240.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0240_singlet_barcodes.txt",header = FALSE)$V1
SNS0691.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0691_singlet_barcodes.txt",header = FALSE)$V1
SNS0786.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/NUC0786_singlet_barcodes.txt",header = FALSE)$V1

all.SCS.singlet.barcodes <- c(SNS0008.singlet.barcodes,
                              SNS010A.singlet.barcodes,
                              SNS015A.singlet.barcodes,
                              SNS0240.singlet.barcodes,
                              SNS0691.singlet.barcodes,
                              SNS0786.singlet.barcodes)

merge$DF.status <- ifelse(rownames(merge@meta.data) %in% all.SCS.singlet.barcodes,"Singlet","Doublet")

merge_singlets <- subset(merge, subset = DF.status == "Singlet")

sample.list <- SplitObject(merge_singlets, split.by = "orig.ident")

sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)
SNS.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
SNS.combined <- IntegrateData(anchorset = SNS.anchors)

DefaultAssay(SNS.combined) <- "integrated"
SNS.combined <- ScaleData(SNS.combined, verbose = FALSE)
SNS.combined <- RunPCA(SNS.combined, npcs = 30, verbose = FALSE)
SNS.combined <- RunUMAP(SNS.combined, reduction = "pca", dims = 1:30)
SNS.combined <- FindNeighbors(SNS.combined, reduction = "pca", dims = 1:30)
SNS.combined <- FindClusters(SNS.combined, resolution = 0.4)

DefaultAssay(SNS.combined) <- "RNA"
SNS.combined <- NormalizeData(SNS.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
SNS.combined <- ScaleData(SNS.combined)

SNS.combined.markers <- FindAllMarkers(SNS.combined, assay = "RNA",logfc.threshold = 0.5,min.pct = 0.25,only.pos = TRUE) %>% filter(p_val_adj<0.05)
write.table(SNS.combined.markers,file = "./Integration/DEG_table/SNS_atlas_DEG_cluster.txt", sep = "\t",quote=FALSE, row.names=FALSE)

new.cluster.ids <- c(
  'CD8', #0
  
  'Hepatocytes', #1
  'Hepatocytes', #2
  'Endothelial', #3
  'Myeloid', #4
  'Hepatocytes', #5
  
  'Cholangiocytes', #6
  'Myeloid', #7
  'Undefine_1', #8
  'CD4', #9
  'Proliferating', #10
  
  'Stellate/Fibroblasts', #11
  'Innate', #12
  'Endothelial', #13
  'Plasma', #14
  'BICC1-hi_hepato', #15
  
  'B', #16
  'Endothelial'#17
)

names(new.cluster.ids) <- levels(SNS.combined)
SNS.combined <- RenameIdents (SNS.combined, new.cluster.ids)
SNS.combined$annotation <- SNS.combined@active.ident

df <- data.frame(SNS_annotation = SNS.combined$annotation,SNS_UMAP1 = SNS.combined@reductions$umap@cell.embeddings[,1], SNS_UMAP2 = SNS.combined@reductions$umap@cell.embeddings[,2])
head(df)
write.csv(df,"./Integration/Table/SNS_atlas_projection_annotation.csv",quote=FALSE)

Idents(SNS.combined) <- "annotation"
levels(SNS.combined) <- c("CD8", "CD4", "Innate", "Myeloid", "B", "Plasma", "Hepatocytes", "BICC1-hi_hepato", "Endothelial", "Cholangiocytes", "Stellate/Fibroblasts", "Proliferating", "Undefine_1")

annotations.markers <- FindAllMarkers(SNS.combined, assay = "RNA",logfc.threshold = 0.5,min.pct = 0.25,only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(annotations.markers,file = "./Integration/DEG_table/SNS_atlas_DEG_annotation.txt", sep = "\t",quote=FALSE, row.names=FALSE)
saveRDS(SNS.combined, file='./Integration/RDS/6_SNS_atlas.rds')