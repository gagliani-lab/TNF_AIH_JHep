## Integration + Annotation Hepatocytes ##

library(dplyr)
library(Seurat)
library(ggplot2)

AIH.combined <- readRDS(file='./Integration/RDS/3_SCS_SNS_atlas.rds')
hepatocytes <- subset(x = AIH.combined, idents = "Hepatocytes")

heaptocytes.SNS <- subset(x = hepatocytes, subset = protocol=="SNS")
DefaultAssay(heaptocytes.SNS) <- "RNA"

sample.list <- SplitObject(heaptocytes.SNS, split.by = "orig.ident")
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)

hepatocytes.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
hepatocytes.combined <- IntegrateData(anchorset = hepatocytes.anchors, dims = 1:30)

DefaultAssay(hepatocytes.combined) <- "integrated"
hepatocytes.combined <- ScaleData(hepatocytes.combined, verbose = FALSE)
hepatocytes.combined <- RunPCA(hepatocytes.combined, npcs = 30, verbose = FALSE)
hepatocytes.combined <- RunUMAP(hepatocytes.combined, reduction = "pca", dims = 1:30)
hepatocytes.combined <- FindNeighbors(hepatocytes.combined, reduction = "pca", dims = 1:30)

hepatocytes.combined <- FindClusters(hepatocytes.combined, resolution = 0.4)

hepatocytes.combined$nCount_CITE <- NULL
hepatocytes.combined$nFeature_CITE <- NULL
hepatocytes.combined$Feature_Count_ratio <- NULL

DefaultAssay(hepatocytes.combined) <- "RNA"
hepatocytes.combined <- NormalizeData(hepatocytes.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
hepatocytes.combined <- ScaleData(hepatocytes.combined) 

hepatocytes.combined.markers <- FindAllMarkers(hepatocytes.combined, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.25, only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(hepatocytes.combined.markers,file = "./Integration/DEG_table/SNS_hepatocytes_DEG_cluster_rev_04.txt", sep = "\t",quote=FALSE, row.names=FALSE)

new.cluster.ids <- c(
  'Midzone', #0
  
  'Portal', #1
  'Central', #2
  'Un-define', #3
  'Inflammatory', #4
  'Proliferation-1', #5
  
  'PTPRC-hi', #6
  'BICC1-hi', #7
  'PTPRB-hi', #8
  'PTPRC-hi', #9
  'Proliferation-2' #10
)

names(new.cluster.ids) <- levels(hepatocytes.combined)
hepatocytes.combined <- RenameIdents (hepatocytes.combined, new.cluster.ids)
hepatocytes.combined$annotation <- hepatocytes.combined@active.ident



levels(hepatocytes.combined)

levels (hepatocytes.combined) <- c("Portal", "Midzone", "Central", "Inflammatory", "Proliferation-1", "Proliferation-2", "BICC1-hi", "PTPRB-hi", "PTPRC-hi", "Un-define")
levels (hepatocytes.combined)

DefaultAssay (hepatocytes.combined)

hepatocytes.combined.annotation.markers <- FindAllMarkers(hepatocytes.combined, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.25, only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(hepatocytes.combined.annotation.markers,file = "./Integration/DEG_table/SNS_hepatocytes_annotation_DEG.txt", sep = "\t",quote=FALSE, row.names=FALSE)

df <- data.frame(annotation = hepatocytes.combined$annotation)
head(df)

write.csv(df,"./Integration/Interactom_annotation/Interactom_hepatocytes_annotation.csv",quote=FALSE)

saveRDS(hepatocytes.combined, file='./Integration/RDS/14_SNS_hepatocytes_figures.rds')