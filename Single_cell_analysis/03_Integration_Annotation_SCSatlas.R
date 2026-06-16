## Integration + Annotation SCS atlas ##

library(dplyr)
library(Seurat)
library(ggplot2)

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


merge <- merge(SCS_05, y = c(SCS_07, 
                             SCS_08, 
                             SCS_09, 
                             SCS_10, 
                             SCS_11, 
                             SCS_12, 
                             SCS_13, 
                             SCS_14, 
                             SCS_15), 
               add.cell.ids = c("AIH05", 
                                "AIH07", 
                                "AIH08", 
                                "AIH09", 
                                "AIH10", 
                                "AIH11", 
                                "AIH12", 
                                "AIH13", 
                                "AIH14", 
                                "AIH15"))


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

all.SCS.singlet.barcodes <- c(AIH05.singlet.barcodes, 
                              AIH07.singlet.barcodes,
                              AIH08.singlet.barcodes,
                              AIH09.singlet.barcodes,
                              AIH10.singlet.barcodes,
                              AIH11.singlet.barcodes,
                              AIH12.singlet.barcodes,
                              AIH13.singlet.barcodes,
                              AIH14.singlet.barcodes,
                              AIH15.singlet.barcodes)

merge$DF.status <- ifelse(rownames(merge@meta.data) %in% all.SCS.singlet.barcodes,"Singlet","Doublet")

merge_singlets <- subset(merge, subset = DF.status == "Singlet")

sample.list <- SplitObject(merge_singlets, split.by = "orig.ident")

sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)

SCS.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)

SCS.combined <- IntegrateData(anchorset = SCS.anchors)

DefaultAssay(SCS.combined) <- "integrated"
SCS.combined <- ScaleData(SCS.combined, verbose = FALSE)
SCS.combined <- RunPCA(SCS.combined, npcs = 30, verbose = FALSE)
SCS.combined <- RunUMAP(SCS.combined, reduction = "pca", dims = 1:30)
SCS.combined <- FindNeighbors(SCS.combined, reduction = "pca", dims = 1:30)
SCS.combined <- FindClusters(SCS.combined, resolution = 0.3)

DefaultAssay(SCS.combined) <- "RNA"
SCS.combined <- NormalizeData(SCS.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
SCS.combined <- ScaleData(SCS.combined)

DefaultAssay(SCS.combined) <- "CITE"
SCS.combined <- NormalizeData(SCS.combined, normalization.method = "CLR", margin = 2, assay = "CITE")
SCS.combined <- ScaleData(SCS.combined)


SCS.combined.markers <- FindAllMarkers(SCS.combined, assay = "RNA",logfc.threshold = 0.5,min.pct = 0.25,only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(SCS.combined.markers,file = "./Integration/DEG_table/SCS_atlas_DEG_cluster.txt", sep = "\t",quote=FALSE, row.names=FALSE)

new.cluster.ids <- c(
  'CD8', #0
  
  'CD4', #1
  'Innate_1', #2
  'Innate_2', #3  
  'Innate_1', #4
  'CD8', #5
  
  'Myeloid', #6
  'MT-hi', #7
  'Innate_like_T', #8
  'B', #9
  'CD4', #10
  
  'Proliferating', #11
  'CD8', #12
  'Innate_like_T', #13
  'Innate_1', #14
  'Plasma_cell' #15
)

names(new.cluster.ids) <- levels(SCS.combined)
SCS.combined <- RenameIdents (SCS.combined, new.cluster.ids)
SCS.combined$annotation <- SCS.combined@active.ident

df <- data.frame(SCS_annotation = SCS.combined$annotation, SCS_UMAP1 = SCS.combined@reductions$umap@cell.embeddings[,1], SCS_UMAP2 = SCS.combined@reductions$umap@cell.embeddings[,2])
write.csv(df,"./Integration/Table/SCS_atlas_projection_annotation.csv",quote=FALSE)

Idents(SCS.combined) <- "annotation"
levels(SCS.combined) <- c("CD8", "CD4", "Myeloid", "Innate_1", "Innate_2", "Innate_like_T","B",
                             "Plasma_cell", "Proliferating", "MT-hi")
annotations.markers <- FindAllMarkers(SCS.combined, assay = "RNA",logfc.threshold = 0.5,min.pct = 0.25,only.pos = TRUE) %>% filter(p_val_adj<0.05)
write.table(annotations.markers,file = "./Integration/DEG_table/SCS_atlas_DEG_annotation.txt", sep = "\t",quote=FALSE, row.names=FALSE)

saveRDS(SCS.combined, file='./Integration/RDS/5_SCS_atlas.rds')
