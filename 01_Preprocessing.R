## Preprocessing of AIH SCS/SNS samples before integration ##

library(dplyr)
library(Seurat)
library(ggplot2)
library(DoubletFinder)



# SCS Liver sample A05 ---------------------------------------------------------

A05 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH5_GEX_ADT/outs/filtered_feature_bc_matrix")

A05 <- CreateA05ratObject(counts = A05$"Gene Expression", project = "A05", 
                          min.cells = 0, 
                          min.features = 0)
A05[["CITE"]] <- CreateAssayObject(counts = A05$"Antibody Capture")
A05[["percent.mt"]] <- PercentageFeatureSet(A05, pattern = "^MT-")
plot1 <- FeatureScatter(A05, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A05, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
A05@meta.data$Feature_Count_ratio <- A05@meta.data$nFeature_RNA/A05@meta.data$nCount_RNA
hist(A05@meta.data$Feature_Count_ratio)
A05 <- subset(A05, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio > 0.2)

saveRDS(A05, file="./Single_sample_preprocess/RDS/prefiltered_AIH05.rds")

A05 <- A05 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A05@var.genes, npc=30)
ElbowPlot(A05)
A05 <- FindNeighbors(A05, dims = 1:15) %>% FindClusters(A05, resolution = 0.5) %>% RunUMAP(A05, dims = 1:15)
DimPlot(A05, reduction = "umap", label=TRUE)


sweep.res.list_liver <- paramSweep_v3(A05, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

homotypic.prop <- modelHomotypic(A05$A05rat_clusters)
nExp_poi <- round(0.007*nrow(A05@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A05 <- doubletFinder_v3(A05, PCs = 1:15, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A05)

Idents(A05) <- "DF.classifications_0.25_0.04_4"
DimPlot(A05, reduction = "umap")

table(A05$DF.classifications_0.25_0.04_4)

A05$barcode <- paste0("AIH05_",rownames(A05@meta.data))
singlet.barcodes <- as.character(A05$barcode[A05$DF.classifications_0.25_0.04_4 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH05_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)

# SCS Liver sample A07 ---------------------------------------------------------

A07 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH7_GEX_ADT/outs/filtered_feature_bc_matrix")

A07 <- CreateSeuratObject(counts = A07$"Gene Expression", project = "A07", 
                          min.cells = 0, 
                          min.features = 0)
A07[["CITE"]] <- CreateAssayObject(counts = A07$"Antibody Capture")
A07[["percent.mt"]] <- PercentageFeatureSet(A07, pattern = "^MT-")

plot1 <- FeatureScatter(A07, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A07, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
A07@meta.data$Feature_Count_ratio <- A07@meta.data$nFeature_RNA/A07@meta.data$nCount_RNA
hist(A07@meta.data$Feature_Count_ratio)

A07 <- subset(A07, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio > 0.2)

saveRDS(A07, file="./Single_sample_preprocess/RDS/prefiltered_AIH07.rds")

A07 <- A07 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A07@var.genes, npc=30)
ElbowPlot(A07)
A07 <- FindNeighbors(A07, dims = 1:15) %>% FindClusters(A07, resolution = 0.5) %>% RunUMAP(A07, dims = 1:15)
DimPlot(A07, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A07, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A07$A07rat_clusters)
nExp_poi <- round(0.016*nrow(A07@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A07 <- doubletFinder_v3(A07, PCs = 1:15, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A07)

Idents(A07) <- "DF.classifications_0.25_0.07_20"
DimPlot(A07, reduction = "umap")

table(A07$DF.classifications_0.25_0.07_20)

A07$barcode <- paste0("AIH07_",rownames(A07@meta.data))
singlet.barcodes <- as.character(A07$barcode[A07$DF.classifications_0.25_0.07_20 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH07_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)

# SCS Liver sample A08 ---------------------------------------------------------

A08 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH8_GEX_ADT/outs/filtered_feature_bc_matrix")

A08 <- CreateSeuratObject(counts = A08$"Gene Expression", project = "A08", 
                          min.cells = 0, 
                          min.features = 0)
A08[["CITE"]] <- CreateAssayObject(counts = A08$"Antibody Capture")
A08[["percent.mt"]] <- PercentageFeatureSet(A08, pattern = "^MT-")

plot1 <- FeatureScatter(A08, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A08, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A08@meta.data$Feature_Count_ratio <- A08@meta.data$nFeature_RNA/A08@meta.data$nCount_RNA
hist(A08@meta.data$Feature_Count_ratio)

A08 <- subset(A08, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)
saveRDS(A08, file="./Single_sample_preprocess/RDS/prefiltered_AIH08.rds")

A08 <- A08 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A08@var.genes, npc=30)
ElbowPlot(A08)

A08 <- FindNeighbors(A08, dims = 1:15) %>% FindClusters(A08, resolution = 0.5) %>% RunUMAP(A08, dims = 1:15)
DimPlot(A08, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A08, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
abline(v=pK_choose[6],lwd=2,col='red',lty=2)
text(pK_choose[6],max(BCmetric),as.character(pK_choose[6]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A08$A08rat_clusters)
nExp_poi <- round(0.023*nrow(A08@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A08 <- doubletFinder_v3(A08, PCs = 1:15, pN = 0.25, pK = 0.21, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A08)

Idents(A08) <- "DF.classifications_0.25_0.21_39"
DimPlot(A08, reduction = "umap")

table(A08$DF.classifications_0.25_0.21_39)

A08$barcode <- paste0("AIH08_",rownames(A08@meta.data))
singlet.barcodes <- as.character(A08$barcode[A08$DF.classifications_0.25_0.21_39 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH08_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A09 ---------------------------------------------------------

A09 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH9_GEX_ADT/filtered_feature_bc_matrix")

A09.adt <- A09$"Antibody Capture"
rownames(A09.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A09.adt))
rownames(A09.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A09.adt))
rownames(A09.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A09.adt))
rownames(A09.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A09.adt))
rownames(A09.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A09.adt))
rownames(A09.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A09.adt))
rownames(A09.adt)<-sub("ADT-TCRab","ADT-TCRAB",rownames(A09.adt))

A09 <- CreateSeuratObject(counts = A09$"Gene Expression", project = "A09", 
                          min.cells = 0, 
                          min.features = 0)
A09[["CITE"]] <- CreateAssayObject(counts = A09.adt)
A09[["percent.mt"]] <- PercentageFeatureSet(A09, pattern = "^MT-")

plot1 <- FeatureScatter(A09, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A09, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A09@meta.data$Feature_Count_ratio <- A09@meta.data$nFeature_RNA/A09@meta.data$nCount_RNA
hist(A09@meta.data$Feature_Count_ratio)

A09 <- subset(A09, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)
saveRDS(A09, file="./Single_sample_preprocess/RDS/prefiltered_AIH09.rds")

A09 <- A09 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A09@var.genes, npc=30)
ElbowPlot(A09)

A09 <- FindNeighbors(A09, dims = 1:15) %>% FindClusters(A09, resolution = 0.5) %>% RunUMAP(A09, dims = 1:15)
DimPlot(A09, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A09, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A09$A09rat_clusters)
nExp_poi <- round(0.023*nrow(A09@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A09 <- doubletFinder_v3(A09, PCs = 1:15, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A09)

Idents(A09) <- "DF.classifications_0.25_0.11_49"
DimPlot(A09, reduction = "umap")

table(A09$DF.classifications_0.25_0.11_49)

A09$barcode <- paste0("AIH09_",rownames(A09@meta.data))
singlet.barcodes <- as.character(A09$barcode[A09$DF.classifications_0.25_0.11_49 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH09_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A10 ---------------------------------------------------------

A10 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH10_GEX_ADT/filtered_feature_bc_matrix")

A10.adt <- A10$"Antibody Capture"

rownames(A10.adt)<-sub("ADT_","ADT-",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A10.adt))
rownames(A10.adt)<-sub("ADT-TCR_A_B","ADT-TCRAB",rownames(A10.adt))

A10 <- CreateSeuratObject(counts = A10$"Gene Expression", project = "A10", 
                          min.cells = 0, 
                          min.features = 0)
A10[["CITE"]] <- CreateAssayObject(counts = A10.adt)

A10[["percent.mt"]] <- PercentageFeatureSet(A10, pattern = "^MT-")

plot1 <- FeatureScatter(A10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A10@meta.data$Feature_Count_ratio <- A10@meta.data$nFeature_RNA/A10@meta.data$nCount_RNA
hist(A10@meta.data$Feature_Count_ratio)

A10 <- subset(A10, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A10, file="./Single_sample_preprocess/RDS/prefiltered_AIH10.rds")

A10 <- A10 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A10@var.genes, npc=30)

ElbowPlot(A10)

A10 <- FindNeighbors(A10, dims = 1:15) %>% FindClusters(A10, resolution = 0.5) %>% RunUMAP(A10, dims = 1:15)
DimPlot(A10, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A10, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
abline(v=pK_choose[6],lwd=2,col='red',lty=2)
text(pK_choose[6],max(BCmetric),as.character(pK_choose[6]),pos = 4,col = "red")
abline(v=pK_choose[7],lwd=2,col='red',lty=2)
text(pK_choose[7],max(BCmetric),as.character(pK_choose[7]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A10$A10rat_clusters)
nExp_poi <- round(0.016*nrow(A10@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A10 <- doubletFinder_v3(A10, PCs = 1:15, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A10)

Idents(A10) <- "DF.classifications_0.25_0.19_20"
DimPlot(A10, reduction = "umap")

table(A10$DF.classifications_0.25_0.19_20)

A10$barcode <- paste0("AIH10_",rownames(A10@meta.data))
singlet.barcodes <- as.character(A10$barcode[A10$DF.classifications_0.25_0.19_20 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH10_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A11 ---------------------------------------------------------

A11 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH11_GEX_ADT/filtered_feature_bc_matrix")

A11.adt <- A11$"Antibody Capture"

rownames(A11.adt)<-sub("ADT_","ADT-",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A11.adt))
rownames(A11.adt)<-sub("ADT-TCR_A_B","ADT-TCRAB",rownames(A11.adt))

A11 <- CreateSeuratObject(counts = A11$"Gene Expression", project = "A11", 
                          min.cells = 0, 
                          min.features = 0)
A11[["CITE"]] <- CreateAssayObject(counts = A11.adt)
A11[["percent.mt"]] <- PercentageFeatureSet(A11, pattern = "^MT-")

plot1 <- FeatureScatter(A11, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A11, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A11@meta.data$Feature_Count_ratio <- A11@meta.data$nFeature_RNA/A11@meta.data$nCount_RNA
hist(A11@meta.data$Feature_Count_ratio)

A11 <- subset(A11, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A11, file="./Single_sample_preprocess/RDS/prefiltered_AIH11.rds")

A11 <- A11 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A11@var.genes, npc=30)
ElbowPlot(A11)
A11 <- FindNeighbors(A11, dims = 1:15) %>% FindClusters(A11, resolution = 0.5) %>% RunUMAP(A11, dims = 1:15)
DimPlot(A11, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A11, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A11$A11rat_clusters)
nExp_poi <- round(0.061*nrow(A11@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A11 <- doubletFinder_v3(A11, PCs = 1:15, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A11)

Idents(A11) <- "DF.classifications_0.25_0.1_400"
DimPlot(A11, reduction = "umap")

table(A11$DF.classifications_0.25_0.1_400)

A11$barcode <- paste0("AIH11_",rownames(A11@meta.data))
singlet.barcodes <- as.character(A11$barcode[A11$DF.classifications_0.25_0.1_400 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH11_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A12 ---------------------------------------------------------

A12 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH12_GEX_ADT/filtered_feature_bc_matrix")

A12.adt <- A12$"Antibody Capture"

rownames(A12.adt)<-sub("ADT_","ADT-",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A12.adt))
rownames(A12.adt)<-sub("ADT-TCR_A_B","ADT-TCRAB",rownames(A12.adt))

A12 <- CreateSeuratObject(counts = A12$"Gene Expression", project = "A12", 
                          min.cells = 0, 
                          min.features = 0)
A12[["CITE"]] <- CreateAssayObject(counts = A12.adt)
A12[["percent.mt"]] <- PercentageFeatureSet(A12, pattern = "^MT-")

plot1 <- FeatureScatter(A12, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A12@meta.data$Feature_Count_ratio <- A12@meta.data$nFeature_RNA/A12@meta.data$nCount_RNA
hist(A12@meta.data$Feature_Count_ratio)

A12 <- subset(A12, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A12, file="./Single_sample_preprocess/RDS/prefiltered_AIH12.rds")

A12 <- A12 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A12@var.genes, npc=30)
ElbowPlot(A12)
A12 <- FindNeighbors(A12, dims = 1:15) %>% FindClusters(A12, resolution = 0.5) %>% RunUMAP(A12, dims = 1:15)
DimPlot(A12, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A12, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A12$A12rat_clusters)
nExp_poi <- round(0.008*nrow(A12@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A12 <- doubletFinder_v3(A12, PCs = 1:15, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A12)

Idents(A12) <- "DF.classifications_0.25_0.05_7"
DimPlot(A12, reduction = "umap", )

table(A12$DF.classifications_0.25_0.05_7)

A12$barcode <- paste0("AIH12_",rownames(A12@meta.data))
singlet.barcodes <- as.character(A12$barcode[A12$DF.classifications_0.25_0.05_7 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH12_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A13 ---------------------------------------------------------

A13 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH13_GEX_ADT/filtered_feature_bc_matrix")

A13.adt <- A13$"Antibody Capture"
rownames(A13.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A13.adt))
rownames(A13.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A13.adt))
rownames(A13.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A13.adt))
rownames(A13.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A13.adt))
rownames(A13.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A13.adt))
rownames(A13.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A13.adt))
rownames(A13.adt)<-sub("ADT-TCRab","ADT-TCRAB",rownames(A13.adt))

A13 <- CreateSeuratObject(counts = A13$"Gene Expression", project = "A13", 
                          min.cells = 0, 
                          min.features = 0)
A13[["CITE"]] <- CreateAssayObject(counts = A13.adt)
A13[["percent.mt"]] <- PercentageFeatureSet(A13, pattern = "^MT-")

plot1 <- FeatureScatter(A13, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A13@meta.data$Feature_Count_ratio <- A13@meta.data$nFeature_RNA/A13@meta.data$nCount_RNA
hist(A13@meta.data$Feature_Count_ratio)

A13 <- subset(A13, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A13, file="./Single_sample_preprocess/RDS/prefiltered_AIH13.rds")

A13 <- A13 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A13@var.genes, npc=30)
ElbowPlot(A13)
A13 <- FindNeighbors(A13, dims = 1:15) %>% FindClusters(A13, resolution = 0.5) %>% RunUMAP(A13, dims = 1:15)
DimPlot(A13, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A13, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A13$A13rat_clusters)
nExp_poi <- round(0.061*nrow(A13@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A13 <- doubletFinder_v3(A13, PCs = 1:15, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A13)

Idents(A13) <- "DF.classifications_0.25_0.08_445"
DimPlot(A13, reduction = "umap")

table(A13$DF.classifications_0.25_0.08_445)

A13$barcode <- paste0("AIH13_",rownames(A13@meta.data))
singlet.barcodes <- as.character(A13$barcode[A13$DF.classifications_0.25_0.08_445 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH13_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A14 ---------------------------------------------------------

A14 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH14_GEX_ADT/filtered_feature_bc_matrix")

A14.adt <- A14$"Antibody Capture"
rownames(A14.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A14.adt))
rownames(A14.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A14.adt))
rownames(A14.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A14.adt))
rownames(A14.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A14.adt))
rownames(A14.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A14.adt))
rownames(A14.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A14.adt))
rownames(A14.adt)<-sub("ADT-TCRab","ADT-TCRAB",rownames(A14.adt))

A14 <- CreateSeuratObject(counts = A14$"Gene Expression", project = "A14", 
                          min.cells = 0, 
                          min.features = 0)
A14[["CITE"]] <- CreateAssayObject(counts = A14.adt)
A14[["percent.mt"]] <- PercentageFeatureSet(A14, pattern = "^MT-")

plot1 <- FeatureScatter(A14, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A14@meta.data$Feature_Count_ratio <- A14@meta.data$nFeature_RNA/A14@meta.data$nCount_RNA
hist(A14@meta.data$Feature_Count_ratio)

A14 <- subset(A14, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A14, file="./Single_sample_preprocess/RDS/prefiltered_AIH14.rds")

A14 <- A14 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A14@var.genes, npc=30)
ElbowPlot(A14)
A14 <- FindNeighbors(A14, dims = 1:15) %>% FindClusters(A14, resolution = 0.5) %>% RunUMAP(A14, dims = 1:15)
DimPlot(A14, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A14, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A14$A14rat_clusters)
nExp_poi <- round(0.054*nrow(A14@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A14 <- doubletFinder_v3(A14, PCs = 1:15, pN = 0.25, pK = 0.23, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A14)

Idents(A14) <- "DF.classifications_0.25_0.23_320"
DimPlot(A14, reduction = "umap")

table(A14$DF.classifications_0.25_0.23_320)

A14$barcode <- paste0("AIH14_",rownames(A14@meta.data))
singlet.barcodes <- as.character(A14$barcode[A14$DF.classifications_0.25_0.23_320 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH14_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS Liver sample A15 ---------------------------------------------------------

A15 <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH15_GEX_ADT/filtered_feature_bc_matrix")

A15.adt <- A15$"Antibody Capture"
rownames(A15.adt)<-sub("ADT_","ADT-",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD197","ADT-CCR7",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD366","ADT-TIM3",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-TCR_A_B","ADT-TCRAB",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD279","ADT-PDCD1",rownames(A15.adt))
rownames(A15.adt)<-sub("ADT-CD152","ADT-CTLA4",rownames(A15.adt))

A15 <- CreateSeuratObject(counts = A15$"Gene Expression", project = "A15", 
                          min.cells = 0, 
                          min.features = 0)
A15[["CITE"]] <- CreateAssayObject(counts = A15.adt)
A15[["percent.mt"]] <- PercentageFeatureSet(A15, pattern = "^MT-")

plot1 <- FeatureScatter(A15, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A15, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A15@meta.data$Feature_Count_ratio <- A15@meta.data$nFeature_RNA/A15@meta.data$nCount_RNA
hist(A15@meta.data$Feature_Count_ratio)

A15 <- subset(A15, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A15, file="./Single_sample_preprocess/RDS/prefiltered_AIH15.rds")

A15 <- A15 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A15@var.genes, npc=30)
ElbowPlot(A15)

A15 <- FindNeighbors(A15, dims = 1:15) %>% FindClusters(A15, resolution = 0.5) %>% RunUMAP(A15, dims = 1:15)
DimPlot(A15, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A15, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A15$A15rat_clusters)
nExp_poi <- round(0.054*nrow(A15@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A15 <- doubletFinder_v3(A15, PCs = 1:15, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A15)

Idents(A15) <- "DF.classifications_0.25_0.11_326"
DimPlot(A15, reduction = "umap")

table(A15$DF.classifications_0.25_0.11_326)

A15$barcode <- paste0("AIH15_",rownames(A15@meta.data))
singlet.barcodes <- as.character(A15$barcode[A15$DF.classifications_0.25_0.11_326 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/AIH15_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SNS Liver sample NUC0008 -----------------------------------------------------

NUC0008.data <- Read10X(data.dir = "./Single_sample_Cellranger/SN-Seq-Cellranger/NucSeq_AIH0008/outs/filtered_feature_bc_matrix")

NUC0008 <- CreateSeuratObject(counts = NUC0008.data, project = "NUC0008", 
                              min.cells = 0, 
                              min.features = 0)
NUC0008[["percent.mt"]] <- PercentageFeatureSet(NUC0008, pattern = "^MT-")

plot1 <- FeatureScatter(NUC0008, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NUC0008, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NUC0008 <- subset(NUC0008, subset = nFeature_RNA > 350 & percent.mt < 1)

saveRDS(NUC0008, file="./Single_sample_preprocess/RDS/prefiltered_NUC0008.rds")

NUC0008 <- NUC0008 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=NUC0008@var.genes, npc=30)
ElbowPlot(NUC0008)
NUC0008 <- FindNeighbors(NUC0008, dims = 1:15) %>% FindClusters(NUC0008, resolution = 0.5) %>% RunUMAP(NUC0008, dims = 1:15)
DimPlot(NUC0008, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(NUC0008, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(NUC0008$NUC0008rat_clusters)
nExp_poi <- round(0.076*nrow(NUC0008@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

NUC0008 <- doubletFinder_v3(NUC0008, PCs = 1:15, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(NUC0008)

Idents(NUC0008) <- "DF.classifications_0.25_0.04_710"
DimPlot(NUC0008, reduction = "umap")

table(NUC0008$DF.classifications_0.25_0.04_710)

NUC0008$barcode <- paste0("NUC0008_",rownames(NUC0008@meta.data))
singlet.barcodes <- as.character(NUC0008$barcode[NUC0008$DF.classifications_0.25_0.04_710 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/NUC0008_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SNS Liver sample NUC010A -----------------------------------------------------

NUC010A.data <- Read10X(data.dir = "./Single_sample_Cellranger/SN-Seq-Cellranger/NucSeq_AIH010A/outs/filtered_feature_bc_matrix")

NUC010A <- CreateSeuratObject(counts = NUC010A.data, project = "NUC010A", 
                              min.cells = 0, 
                              min.features = 0)
NUC010A[["percent.mt"]] <- PercentageFeatureSet(NUC010A, pattern = "^MT-")

plot1 <- FeatureScatter(NUC010A, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NUC010A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NUC010A <- subset(NUC010A, subset = nFeature_RNA > 350 & percent.mt < 1 )

saveRDS(NUC010A, file="./Single_sample_preprocess/RDS/prefiltered_NUC010A.rds")

NUC010A <- NUC010A %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=NUC010A@var.genes, npc=30)
ElbowPlot(NUC010A)
NUC010A <- FindNeighbors(NUC010A, dims = 1:15) %>% FindClusters(NUC010A, resolution = 0.5) %>% RunUMAP(NUC010A, dims = 1:15)
DimPlot(NUC010A, reduction = "umap", label=TRUE)

head(NUC010A)

sweep.res.list_liver <- paramSweep_v3(NUC010A, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(NUC010A$NUC010Arat_clusters)
nExp_poi <- round(0.023*nrow(NUC010A@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

NUC010A <- doubletFinder_v3(NUC010A, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(NUC010A)

Idents(NUC010A) <- "DF.classifications_0.25_0.02_58"
DimPlot(NUC010A, reduction = "umap")

table(NUC010A$DF.classifications_0.25_0.02_58)

NUC010A$barcode <- paste0("NUC010A_",rownames(NUC010A@meta.data))
singlet.barcodes <- as.character(NUC010A$barcode[NUC010A$DF.classifications_0.25_0.02_58 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/NUC010A_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SNS Liver sample NUC015A -----------------------------------------------------

NUC015A.data <- Read10X(data.dir = "./Single_sample_Cellranger/SN-Seq-Cellranger/NucSeq_AIH015A/outs/filtered_feature_bc_matrix")

NUC015A <- CreateSeuratObject(counts = NUC015A.data, project = "NUC015A", 
                              min.cells = 0, 
                              min.features = 0)
NUC015A[["percent.mt"]] <- PercentageFeatureSet(NUC015A, pattern = "^MT-")

plot1 <- FeatureScatter(NUC015A, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NUC015A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NUC015A <- subset(NUC015A, subset = nFeature_RNA > 350 & percent.mt < 1)

saveRDS(NUC015A, file="./Single_sample_preprocess/RDS/prefiltered_NUC015A.rds")

NUC015A <- NUC015A %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=NUC015A@var.genes, npc=30)
ElbowPlot(NUC015A)
NUC015A <- FindNeighbors(NUC015A, dims = 1:15) %>% FindClusters(NUC015A, resolution = 0.5) %>% RunUMAP(NUC015A, dims = 1:15)
DimPlot(NUC015A, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(NUC015A, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
abline(v=pK_choose[6],lwd=2,col='red',lty=2)
text(pK_choose[6],max(BCmetric),as.character(pK_choose[6]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(NUC015A$NUC015Arat_clusters)
nExp_poi <- round(0.046*nrow(NUC015A@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

NUC015A <- doubletFinder_v3(NUC015A, PCs = 1:15, pN = 0.25, pK = 0.14, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(NUC015A)

Idents(NUC015A) <- "DF.classifications_0.25_0.14_260"
DimPlot(NUC015A, reduction = "umap")

table(NUC015A$DF.classifications_0.25_0.14_260)

NUC015A$barcode <- paste0("NUC015A_",rownames(NUC015A@meta.data))
singlet.barcodes <- as.character(NUC015A$barcode[NUC015A$DF.classifications_0.25_0.14_260 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/NUC015A_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SNS Liver sample NUC0240 -----------------------------------------------------

NUC0240.data <- Read10X(data.dir = "./Single_sample_Cellranger/SN-Seq-Cellranger/NucSeq_AIH0240/outs/filtered_feature_bc_matrix")

NUC0240 <- CreateSeuratObject(counts = NUC0240.data, project = "NUC0240", 
                              min.cells = 0, 
                              min.features = 0)
NUC0240[["percent.mt"]] <- PercentageFeatureSet(NUC0240, pattern = "^MT-")

plot1 <- FeatureScatter(NUC0240, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NUC0240, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NUC0240 <- subset(NUC0240, subset = nFeature_RNA > 350 & percent.mt < 1)

saveRDS(NUC0240, file="./Single_sample_preprocess/RDS/prefiltered_NUC0240.rds")

NUC0240 <- NUC0240 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=NUC0240@var.genes, npc=30)
ElbowPlot(NUC0240)
NUC0240 <- FindNeighbors(NUC0240, dims = 1:15) %>% FindClusters(NUC0240, resolution = 0.5) %>% RunUMAP(NUC0240, dims = 1:15)
DimPlot(NUC0240, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(NUC0240, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(NUC0240$NUC0240rat_clusters)
nExp_poi <- round(0.092*nrow(NUC0240@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

NUC0240 <- doubletFinder_v3(NUC0240, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(NUC0240)

Idents(NUC0240) <- "DF.classifications_0.25_0.03_1020"
DimPlot(NUC0240, reduction = "umap")

table(NUC0240$DF.classifications_0.25_0.03_1020)

NUC0240$barcode <- paste0("NUC0240_",rownames(NUC0240@meta.data))
singlet.barcodes <- as.character(NUC0240$barcode[NUC0240$DF.classifications_0.25_0.03_1020 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/NUC0240_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SNS Liver sample NUC0691 -----------------------------------------------------

NUC0691.data <- Read10X(data.dir = "./Single_sample_Cellranger/SN-Seq-Cellranger/NucSeq_AIH0691/outs/filtered_feature_bc_matrix")

NUC0691 <- CreateSeuratObject(counts = NUC0691.data, project = "NUC0691", 
                              min.cells = 0, 
                              min.features = 0)
NUC0691[["percent.mt"]] <- PercentageFeatureSet(NUC0691, pattern = "^MT-")

plot1 <- FeatureScatter(NUC0691, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NUC0691, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NUC0691 <- subset(NUC0691, subset = nFeature_RNA > 350 & percent.mt < 1)

saveRDS(NUC0691, file="./Single_sample_preprocess/RDS/prefiltered_NUC0691.rds")

NUC0691 <- NUC0691 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=NUC0691@var.genes, npc=30)
ElbowPlot(NUC0691)
NUC0691 <- FindNeighbors(NUC0691, dims = 1:15) %>% FindClusters(NUC0691, resolution = 0.5) %>% RunUMAP(NUC0691, dims = 1:15)
DimPlot(NUC0691, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(NUC0691, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
abline(v=pK_choose[6],lwd=2,col='red',lty=2)
text(pK_choose[6],max(BCmetric),as.character(pK_choose[6]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(NUC0691$NUC0691rat_clusters)
nExp_poi <- round(0.069*nrow(NUC0691@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

NUC0691 <- doubletFinder_v3(NUC0691, PCs = 1:15, pN = 0.25, pK = 0.2, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(NUC0691)

Idents(NUC0691) <- "DF.classifications_0.25_0.2_615"
DimPlot(NUC0691, reduction = "umap")

table(NUC0691$DF.classifications_0.25_0.2_615)

NUC0691$barcode <- paste0("NUC0691_",rownames(NUC0691@meta.data))
singlet.barcodes <- as.character(NUC0691$barcode[NUC0691$DF.classifications_0.25_0.2_615 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/NUC0691_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SNS Liver sample NUC0786 -----------------------------------------------------

NUC0786.data <- Read10X(data.dir = "./Single_sample_Cellranger/SN-Seq-Cellranger/NucSeq_AIH0786/outs/filtered_feature_bc_matrix")

NUC0786 <- CreateSeuratObject(counts = NUC0786.data, project = "NUC0786", 
                              min.cells = 0, 
                              min.features = 0)
NUC0786[["percent.mt"]] <- PercentageFeatureSet(NUC0786, pattern = "^MT-")

plot1 <- FeatureScatter(NUC0786, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NUC0786, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

NUC0786 <- subset(NUC0786, subset = nFeature_RNA > 350 & percent.mt < 1)

saveRDS(NUC0786, file="./Single_sample_preprocess/RDS/prefiltered_NUC0786.rds")

NUC0786 <- NUC0786 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=NUC0786@var.genes, npc=30)
ElbowPlot(NUC0786)
NUC0786 <- FindNeighbors(NUC0786, dims = 1:15) %>% FindClusters(NUC0786, resolution = 0.5) %>% RunUMAP(NUC0786, dims = 1:15)
DimPlot(NUC0786, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(NUC0786, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
abline(v=pK_choose[6],lwd=2,col='red',lty=2)
text(pK_choose[6],max(BCmetric),as.character(pK_choose[6]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(NUC0786$NUC0786rat_clusters)
nExp_poi <- round(0.061*nrow(NUC0786@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

NUC0786 <- doubletFinder_v3(NUC0786, PCs = 1:15, pN = 0.25, pK = 0.27, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(NUC0786)

Idents(NUC0786) <- "DF.classifications_0.25_0.27_480"
DimPlot(NUC0786, reduction = "umap")

table(NUC0786$DF.classifications_0.25_0.27_480)

NUC0786$barcode <- paste0("NUC0786_",rownames(NUC0786@meta.data))
singlet.barcodes <- as.character(NUC0786$barcode[NUC0786$DF.classifications_0.25_0.27_480 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/NUC0786_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS PBMC sample A09 ----------------------------------------------------------

A09_PBMC <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH9_PBMC_GEX_ADT/filtered_feature_bc_matrix")

A09_PBMC.adt <- A09_PBMC$"Antibody Capture"

rownames(A09_PBMC.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A09_PBMC.adt))
rownames(A09_PBMC.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A09_PBMC.adt))
rownames(A09_PBMC.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A09_PBMC.adt))
rownames(A09_PBMC.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A09_PBMC.adt))
rownames(A09_PBMC.adt)<-sub("ADT-TCR-alpha-beta","ADT-TCRAB",rownames(A09_PBMC.adt))
rownames(A09_PBMC.adt)<-sub("ADT-CD279","ADT-PDCD1",rownames(A09_PBMC.adt))
rownames(A09_PBMC.adt)<-sub("ADT-CD152","ADT-CTLA4",rownames(A09_PBMC.adt))

A09_PBMC <- CreateSeuratObject(counts = A09_PBMC$"Gene Expression", project = "PBMC-A09", 
                               min.cells = 0, 
                               min.features = 0)
A09_PBMC[["CITE"]] <- CreateAssayObject(counts = A09_PBMC.adt)
A09_PBMC[["percent.mt"]] <- PercentageFeatureSet(A09_PBMC, pattern = "^MT-")

plot1 <- FeatureScatter(A09_PBMC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A09_PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A09_PBMC@meta.data$Feature_Count_ratio <- A09_PBMC@meta.data$nFeature_RNA/A09_PBMC@meta.data$nCount_RNA
hist(A09_PBMC@meta.data$Feature_Count_ratio)

A09_PBMC <- subset(A09_PBMC, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A09_PBMC, file="./Single_sample_preprocess/RDS/prefiltered_PBMC_AIH09.rds")

A09_PBMC <- A09_PBMC %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A09_PBMC@var.genes, npc=30)
ElbowPlot(A09_PBMC)
A09_PBMC <- FindNeighbors(A09_PBMC, dims = 1:15) %>% FindClusters(A09_PBMC, resolution = 0.5) %>% RunUMAP(A09_PBMC, dims = 1:15)
DimPlot(A09_PBMC, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A09_PBMC, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A09_PBMC$A09_PBMCrat_clusters)
nExp_poi <- round(0.085*nrow(A09_PBMC@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A09_PBMC <- doubletFinder_v3(A09_PBMC, PCs = 1:15, pN = 0.25, pK = 0.27, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A09_PBMC)

Idents(A09_PBMC) <- "DF.classifications_0.25_0.27_949"
DimPlot(A09_PBMC, reduction = "umap")

table(A09_PBMC$DF.classifications_0.25_0.27_949)

A09_PBMC$barcode <- paste0("PBMC_AIH09_",rownames(A09_PBMC@meta.data))
singlet.barcodes <- as.character(A09_PBMC$barcode[A09_PBMC$DF.classifications_0.25_0.27_949 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/PBMC_AIH09_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS PBMC sample A13 ----------------------------------------------------------

A13_PBMC <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH13_PBMC_GEX_ADT/filtered_feature_bc_matrix")

A13_PBMC.adt <- A13_PBMC$"Antibody Capture"
rownames(A13_PBMC.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A13_PBMC.adt))
rownames(A13_PBMC.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A13_PBMC.adt))
rownames(A13_PBMC.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A13_PBMC.adt))
rownames(A13_PBMC.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A13_PBMC.adt))
rownames(A13_PBMC.adt)<-sub("ADT-TCR-alpha-beta","ADT-TCRAB",rownames(A13_PBMC.adt))
rownames(A13_PBMC.adt)<-sub("ADT-CD279","ADT-PDCD1",rownames(A13_PBMC.adt))
rownames(A13_PBMC.adt)<-sub("ADT-CD152","ADT-CTLA4",rownames(A13_PBMC.adt))

A13_PBMC <- CreateSeuratObject(counts = A13_PBMC$"Gene Expression", project = "PBMC-A13", 
                               min.cells = 0, 
                               min.features = 0)
A13_PBMC[["CITE"]] <- CreateAssayObject(counts = A13_PBMC.adt)
A13_PBMC[["percent.mt"]] <- PercentageFeatureSet(A13_PBMC, pattern = "^MT-")

plot1 <- FeatureScatter(A13_PBMC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A13_PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A13_PBMC@meta.data$Feature_Count_ratio <- A13_PBMC@meta.data$nFeature_RNA/A13_PBMC@meta.data$nCount_RNA
hist(A13_PBMC@meta.data$Feature_Count_ratio)

A13_PBMC <- subset(A13_PBMC, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A13_PBMC, file="./Single_sample_preprocess/RDS/prefiltered_PBMC_AIH13.rds")

A13_PBMC <- A13_PBMC %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A13_PBMC@var.genes, npc=30)
ElbowPlot(A13_PBMC)
A13_PBMC <- FindNeighbors(A13_PBMC, dims = 1:15) %>% FindClusters(A13_PBMC, resolution = 0.5) %>% RunUMAP(A13_PBMC, dims = 1:15)
DimPlot(A13_PBMC, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A13_PBMC, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A13_PBMC$A13_PBMCrat_clusters)
nExp_poi <- round(0.076*nrow(A13_PBMC@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A13_PBMC <- doubletFinder_v3(A13_PBMC, PCs = 1:15, pN = 0.25, pK = 0.21, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A13_PBMC)

Idents(A13_PBMC) <- "DF.classifications_0.25_0.21_661"
DimPlot(A13_PBMC, reduction = "umap")

table(A13_PBMC$DF.classifications_0.25_0.21_661)

A13_PBMC$barcode <- paste0("PBMC_AIH13_",rownames(A13_PBMC@meta.data))
singlet.barcodes <- as.character(A13_PBMC$barcode[A13_PBMC$DF.classifications_0.25_0.21_661 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/PBMC_AIH13_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)


# SCS PBMC sample A15 ----------------------------------------------------------

A15_PBMC <- Read10X(data.dir = "./Single_sample_Cellranger/SC-Seq-Cellranger/AIH15_PBMC_GEX_ADT/filtered_feature_bc_matrix")

A15_PBMC.adt <- A15_PBMC$"Antibody Capture"
rownames(A15_PBMC.adt)<-sub("ADT-CD196","ADT-CCR6",rownames(A15_PBMC.adt))
rownames(A15_PBMC.adt)<-sub("ADT-CD223","ADT-LAG3",rownames(A15_PBMC.adt))
rownames(A15_PBMC.adt)<-sub("ADT-CD183","ADT-CXCR3",rownames(A15_PBMC.adt))
rownames(A15_PBMC.adt)<-sub("ADT-CD127","ADT-IL7RA",rownames(A15_PBMC.adt))
rownames(A15_PBMC.adt)<-sub("ADT-TCR-alpha-beta","ADT-TCRAB",rownames(A15_PBMC.adt))
rownames(A15_PBMC.adt)<-sub("ADT-CD279","ADT-PDCD1",rownames(A15_PBMC.adt))
rownames(A15_PBMC.adt)<-sub("ADT-CD152","ADT-CTLA4",rownames(A15_PBMC.adt))

A15_PBMC <- CreateSeuratObject(counts = A15_PBMC$"Gene Expression", project = "PBMC-A15", 
                               min.cells = 0, 
                               min.features = 0)
A15_PBMC[["CITE"]] <- CreateAssayObject(counts = A15_PBMC.adt)
A15_PBMC[["percent.mt"]] <- PercentageFeatureSet(A15_PBMC, pattern = "^MT-")

plot1 <- FeatureScatter(A15_PBMC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(A15_PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

A15_PBMC@meta.data$Feature_Count_ratio <- A15_PBMC@meta.data$nFeature_RNA/A15_PBMC@meta.data$nCount_RNA
hist(A15_PBMC@meta.data$Feature_Count_ratio)

A15_PBMC <- subset(A15_PBMC, subset = nFeature_RNA > 350 & percent.mt < 5 & Feature_Count_ratio  > 0.2)

saveRDS(A15_PBMC, file="./Single_sample_preprocess/RDS/prefiltered_PBMC_AIH15.rds")

A15_PBMC <- A15_PBMC %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA (pc.genes=A15_PBMC@var.genes, npc=30)
ElbowPlot(A15_PBMC)
A15_PBMC <- FindNeighbors(A15_PBMC, dims = 1:15) %>% FindClusters(A15_PBMC, resolution = 0.5) %>% RunUMAP(A15_PBMC, dims = 1:15)
DimPlot(A15_PBMC, reduction = "umap", label=TRUE)

sweep.res.list_liver <- paramSweep_v3(A15_PBMC, PCs = 1:15, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_liver)

pK=as.numeric(as.character(bcmvn_meta$pK))
BCmetric=bcmvn_meta$BCmetric
# pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pK_choose = pK[which(diff(sign(diff(BCmetric)))==-2)+1]
pK_choose

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose[1],lwd=2,col='red',lty=2)
text(pK_choose[1],max(BCmetric),as.character(pK_choose[1]),pos = 4,col = "red")
abline(v=pK_choose[2],lwd=2,col='red',lty=2)
text(pK_choose[2],max(BCmetric),as.character(pK_choose[2]),pos = 4,col = "red")
abline(v=pK_choose[3],lwd=2,col='red',lty=2)
text(pK_choose[3],max(BCmetric),as.character(pK_choose[3]),pos = 4,col = "red")
abline(v=pK_choose[4],lwd=2,col='red',lty=2)
text(pK_choose[4],max(BCmetric),as.character(pK_choose[4]),pos = 4,col = "red")
abline(v=pK_choose[5],lwd=2,col='red',lty=2)
text(pK_choose[5],max(BCmetric),as.character(pK_choose[5]),pos = 4,col = "red")
abline(v=pK_choose[6],lwd=2,col='red',lty=2)
text(pK_choose[6],max(BCmetric),as.character(pK_choose[6]),pos = 4,col = "red")
title("The BCmvn distributions")

homotypic.prop <- modelHomotypic(A15_PBMC$A15_PBMCrat_clusters)
nExp_poi <- round(0.061*nrow(A15_PBMC@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

A15_PBMC <- doubletFinder_v3(A15_PBMC, PCs = 1:15, pN = 0.25, pK = 0.12, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(A15_PBMC)

Idents(A15_PBMC) <- "DF.classifications_0.25_0.12_464"
DimPlot(A15_PBMC, reduction = "umap")

table(A15_PBMC$DF.classifications_0.25_0.12_464)

A15_PBMC$barcode <- paste0("PBMC_AIH15_",rownames(A15_PBMC@meta.data))
singlet.barcodes <- as.character(A15_PBMC$barcode[A15_PBMC$DF.classifications_0.25_0.12_464 == "Singlet"])
head(singlet.barcodes) 
length(singlet.barcodes)

write.table(data.frame(barcodes = singlet.barcodes),file="./Single_sample_preprocess/Doublet_finder_table/PBMC_AIH15_singlet_barcodes.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)
