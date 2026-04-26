## Integration + Annotation CD4+ T cells ##

library(dplyr)
library(Seurat)
library(ggplot2)

AIH.combined <- readRDS(file='./Integration/RDS/3_SCS_SNS_atlas.rds')

T_cell <- subset(x = AIH.combined, idents = c("CD4", "CD8"))

sample.list <- SplitObject(T_cell, split.by = "orig.ident")
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)
Tcells.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)

Tcells.combined <- IntegrateData(anchorset = Tcells.anchors, dims = 1:30)

DefaultAssay(Tcells.combined) <- "integrated"
Tcells.combined <- ScaleData(Tcells.combined, verbose = FALSE)
Tcells.combined <- RunPCA(Tcells.combined, npcs = 30, verbose = FALSE)
Tcells.combined <- RunUMAP(Tcells.combined, reduction = "pca", dims = 1:30)
Tcells.combined <- FindNeighbors(Tcells.combined, reduction = "pca", dims = 1:30)
Tcells.combined <- FindClusters(Tcells.combined, resolution = 0.5)

DefaultAssay(Tcells.combined) <- "RNA"
Tcells.combined <- NormalizeData(Tcells.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
Tcells.combined <- ScaleData(Tcells.combined) 

DefaultAssay(Tcells.combined) <- "CITE"
Tcells.combined <- NormalizeData(Tcells.combined, normalization.method = "CLR", margin = 2, assay = "CITE")
Tcells.combined <- ScaleData(Tcells.combined) 

DefaultAssay(Tcells.combined) <- "RNA"
CD8 <- subset(x = Tcells.combined, subset = integrated_snn_res.0.5 %in% c(1,4,8,9,10) & CD8A == 0, invert=TRUE)

CD8$seurat_clusters <- NULL
CD8$annotation <- NULL
CD8$integrated_snn_res.0.5 <- NULL

Idents(CD8) <- "orig.ident"
DefaultAssay(CD8) <- "RNA"

sample.list <- SplitObject(CD8, split.by = "orig.ident")

sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)
CD8.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
CD8.combined <- IntegrateData(anchorset = CD8.anchors, dims = 1:30)

DefaultAssay(CD8.combined) <- "integrated"
CD8.combined <- ScaleData(CD8.combined, verbose = FALSE)
CD8.combined <- RunPCA(CD8.combined, npcs = 30, verbose = FALSE)
CD8.combined <- RunUMAP(CD8.combined, reduction = "pca", dims = 1:30)
CD8.combined <- FindNeighbors(CD8.combined, reduction = "pca", dims = 1:30)
CD8.combined <- FindClusters(CD8.combined, resolution = 0.6)

DefaultAssay(CD8.combined) <- "RNA"
CD8.combined <- NormalizeData(CD8.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
CD8.combined <- ScaleData(CD8.combined) 

DefaultAssay(CD8.combined) <- "CITE"
CD8.combined <- NormalizeData(CD8.combined, normalization.method = "CLR", margin = 2, assay = "CITE")
CD8.combined <- ScaleData(CD8.combined) 

DefaultAssay (CD8.combined) <- "RNA"

CD8.combined.markers <- FindAllMarkers(CD8.combined, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(CD8.combined.markers,file = "./Integration/DEG_table/SCS_SNS_CD8_DEG_cluster.txt", sep = "\t",quote=FALSE, row.names=FALSE)

new.cluster.ids <- c(
  'TEM-d', #0
  
  'TEM-a', #1
  'TEM-e', #2
  'TEM-g', #3
  'TCM-b', #4
  'TCM-a', #5
  
  'TEM-f', #6
  'MT-high', #7
  'TEM-c', #8
  'Naive', #9
  'TEM-b', #10
  
  'CD4' #11
)

names(new.cluster.ids) <- levels(CD8.combined)
CD8.combined <- RenameIdents (CD8.combined, new.cluster.ids)
CD8.combined$annotation <- CD8.combined@active.ident

saveRDS(CD8.combined, file='./Integration/RDS/13_SCS_SNS_CD8_figures.rds')