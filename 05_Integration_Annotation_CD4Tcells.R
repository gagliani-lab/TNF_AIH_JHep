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
CD4 <- subset(x = Tcells.combined, subset = integrated_snn_res.0.5 %in% c(1,4,8,9,10) & CD8A == 0)

CD4$annotation <- NULL
CD4$integrated_snn_res.0.5 <- NULL
CD4$seurat_clusters <- NULL
head(CD4)

saveRDS(CD4, file='./Integration/RDS/10_SCS_SNS_CD4.rds')

# Remove two samples with very few CD4 T cells
CD4_high <- subset(x = CD4, idents = c("A05", "NUC010A"), invert = TRUE)
DefaultAssay(CD4_high) <- "RNA"

sample.list <- SplitObject(CD4_high, split.by = "orig.ident")

sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)
CD4.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
CD4.combined <- IntegrateData(anchorset = CD4.anchors, dims = 1:30)

DefaultAssay(CD4.combined) <- "integrated"
CD4.combined <- ScaleData(CD4.combined, verbose = FALSE)
CD4.combined <- RunPCA(CD4.combined, npcs = 30, verbose = FALSE)
CD4.combined <- RunUMAP(CD4.combined, reduction = "pca", dims = 1:30)
CD4.combined <- FindNeighbors(CD4.combined, reduction = "pca", dims = 1:30)
CD4.combined <- FindClusters(CD4.combined, resolution = 0.7)

Idents(CD4.combined) <- "seurat_clusters"
CD4.combined <- FindSubCluster(CD4.combined, 
                                  "2", 
                                  "integrated_snn", 
                                  resolution = 0.2, 
                                  algorithm = 1)

DefaultAssay(CD4.combined) <- "RNA"
CD4.combined <- NormalizeData(CD4.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
CD4.combined <- ScaleData(CD4.combined) 

DefaultAssay(CD4.combined) <- "CITE"
CD4.combined <- NormalizeData(CD4.combined, normalization.method = "CLR", margin = 2, assay = "CITE")
CD4.combined <- ScaleData(CD4.combined) 

Idents(CD4.combined) <- "sub.cluster"
CD4.combined.markers <- FindAllMarkers(CD4.combined, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE) %>% filter(p_val_adj<0.05)

write.table(CD4.combined.markers,file = "./Integration/DEG_table/SCS_SNS_CD4_DEG_cluster.txt", sep = "\t",quote=FALSE, row.names=FALSE)

new.cluster.ids <- c(
  'TCM', #0
  
  'TEM-b', #2_0
  'Naive', #3
  'TREG', #4
  'TRM1', #1
  'TREG-naive', #7
  
  'TFH', #5
  'TEM-a', #6
  'TR1' #2_1
)

names(new.cluster.ids) <- levels(CD4.combined)
CD4.combined <- RenameIdents (CD4.combined, new.cluster.ids)
CD4.combined$annotation <- CD4.combined@active.ident

## Export barcodes of TRM1 cells for overlay on PBMC/Liver T cells
TRM1.barcodes <- as.character(rownames(CD4.combined@meta.data[CD4.combined$annotation== "TRM1"&(CD4.combined$orig.ident %in% c("A09", "A13", "A15"))]))

write.table(data.frame(barcodes = TRM1.barcodes),file="./Integration/Table/TRM1_barcode.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)

## Add TCR annotation

s7 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH7_TCR/outs/filtered_contig_annotations.csv")
s8 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH8_TCR/outs/filtered_contig_annotations.csv")
s9 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH9_TCR/filtered_contig_annotations.csv")
s10 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH10_TCR/filtered_contig_annotations.csv")
s11 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH11_TCR/filtered_contig_annotations.csv")
s12 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH12_TCR/filtered_contig_annotations.csv")
s13 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH13_TCR/filtered_contig_annotations.csv")
s14 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH14_TCR/filtered_contig_annotations.csv")
s15 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH15_TCR/filtered_contig_annotations.csv")

ob.list <- list ()
samplelist <- c("AIH07", "AIH08", "AIH09", "AIH10", "AIH11", "AIH12", "AIH13", "AIH14", "AIH15")
ob.list[[1]] <- s7
ob.list[[2]] <- s8
ob.list[[3]] <- s9
ob.list[[4]] <- s10
ob.list[[5]] <- s11
ob.list[[6]] <- s12
ob.list[[7]] <- s13
ob.list[[8]] <- s14
ob.list[[9]] <- s15

combined <- combineTCR(ob.list,
                       sample=samplelist,
                       filterMulti=FALSE)

CD4.combined <- combineExpression(combined, 
                         CD4.combined,
                         cloneCall = "gene+nt",
                         proportion = F,
                         cloneType = c(None = 0, Single =1, Small = 3, Medium = 5, Large =20),
                         filterNA= T)



saveRDS(CD4.combined, file='./Integration/RDS/12_SCS_SNS_CD4_figures.rds')


