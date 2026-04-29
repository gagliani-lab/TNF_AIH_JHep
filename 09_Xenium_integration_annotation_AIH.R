## Integration + Annotation Xenium AIH data ##

library(Seurat) # packageVersion("Seurat") # 5.0.3
library(tidyverse) # packageVersion("tidyverse") # 2.0.0
library(harmony) # packageVersion("harmony") # 1.2.0
library(patchwork) # packageVersion("patchwork") # 1.2.0
library(openxlsx) # packageVersion("openxlsx") # 4.2.5.2

## Insert directory here which contains the following folders:
## rawData #(contains Xenium raw data)
## results
## doc
## SeuratObjects
# working.dir <- "path/to/working/directory"

## 1. Import AIH data ----------------------------------------------------------

### 1.1 AIH Pat 1

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_1__20240620__110232/")
AIH.pat1 <- LoadXenium(path, fov = "fov")
AIH.pat1 <- subset(AIH.pat1, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat1),]
csv.cells <- csv.cells[match(colnames(AIH.pat1),csv.cells$cell_id),]
all(colnames(AIH.pat1)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat1 <- AddMetaData(AIH.pat1,csv.cells)
AIH.pat1 <- RenameCells(AIH.pat1,add.cell.id = "AIH.pat1")
AIH.pat1$sample <- "AIH.pat1"


### 1.2 AIH Pat 2 Block 1

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_2__20240620__110232/")
AIH.pat2.1 <- LoadXenium(path, fov = "fov")
AIH.pat2.1 <- subset(AIH.pat2.1, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat2.1),]
csv.cells <- csv.cells[match(colnames(AIH.pat2.1),csv.cells$cell_id),]
all(colnames(AIH.pat2.1)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat2.1 <- AddMetaData(AIH.pat2.1,csv.cells)
AIH.pat2.1 <- RenameCells(AIH.pat2.1,add.cell.id = "AIH.pat2.1")
AIH.pat2.1$sample <- "AIH.pat2.1"


### 1.3 AIH Pat 2 Block 2

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_3__20240620__110232/")
AIH.pat2.2 <- LoadXenium(path, fov = "fov")
AIH.pat2.2 <- subset(AIH.pat2.2, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat2.2),]
csv.cells <- csv.cells[match(colnames(AIH.pat2.2),csv.cells$cell_id),]
all(colnames(AIH.pat2.2)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat2.2 <- AddMetaData(AIH.pat2.2,csv.cells)
AIH.pat2.2 <- RenameCells(AIH.pat2.2,add.cell.id = "AIH.pat2.2")
AIH.pat2.2$sample <- "AIH.pat2.2"


### 1.4 AIH Pat 3

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_4__20240620__110232/")
AIH.pat3 <- LoadXenium(path, fov = "fov")
AIH.pat3 <- subset(AIH.pat3, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat3),]
csv.cells <- csv.cells[match(colnames(AIH.pat3),csv.cells$cell_id),]
all(colnames(AIH.pat3)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat3 <- AddMetaData(AIH.pat3,csv.cells)
AIH.pat3 <- RenameCells(AIH.pat3,add.cell.id = "AIH.pat3")
AIH.pat3$sample <- "AIH.pat3"


### 1.5 AIH Pat 4

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_1__20240620__110232/")
AIH.pat4 <- LoadXenium(path, fov = "fov")
AIH.pat4 <- subset(AIH.pat4, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat4),]
csv.cells <- csv.cells[match(colnames(AIH.pat4),csv.cells$cell_id),]
all(colnames(AIH.pat4)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat4 <- AddMetaData(AIH.pat4,csv.cells)
AIH.pat4 <- RenameCells(AIH.pat4,add.cell.id = "AIH.pat4")
AIH.pat4$sample <- "AIH.pat4"


### 1.6 AIH Pat 5

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_2__20240620__110232/")
AIH.pat5 <- LoadXenium(path, fov = "fov")
AIH.pat5 <- subset(AIH.pat5, subset = nCount_Xenium > 0)

csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat5),]
csv.cells <- csv.cells[match(colnames(AIH.pat5),csv.cells$cell_id),]
all(colnames(AIH.pat5)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat5 <- AddMetaData(AIH.pat5,csv.cells)
AIH.pat5 <- RenameCells(AIH.pat5,add.cell.id = "AIH.pat5")
AIH.pat5$sample <- "AIH.pat5"


### 1.7 AIH Pat 6
## Cut this so only the two pieces at the top are included (the bottom piece is a control biopsy)

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_3__20240620__110232/")
AIH.pat6 <- LoadXenium(path, fov = "fov")
AIH.pat6 <- subset(AIH.pat6, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat6),]
csv.cells <- csv.cells[match(colnames(AIH.pat6),csv.cells$cell_id),]
all(colnames(AIH.pat6)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat6 <- AddMetaData(AIH.pat6,csv.cells)
AIH.pat6 <- subset(AIH.pat6, subset = y_centroid<4500)
AIH.pat6 <- RenameCells(AIH.pat6,add.cell.id = "AIH.pat6")
AIH.pat6$sample <- "AIH.pat6"



## 2. Integrate AIH data -------------------------------------------------------

AIH.combined.integrated <- merge(AIH.pat1,y = c(AIH.pat2.1,AIH.pat2.2,AIH.pat3,AIH.pat4,AIH.pat5,AIH.pat6))

AIH.combined.integrated <- SCTransform(AIH.combined.integrated, assay = "Xenium")
AIH.combined.integrated <- RunPCA(AIH.combined.integrated, npcs = 30, features = rownames(AIH.combined.integrated))
AIH.combined.integrated <- AIH.combined.integrated %>% RunHarmony("sample",plot_convergence = TRUE,max_iter = 15,early_stop=FALSE)
AIH.combined.integrated <- RunUMAP(AIH.combined.integrated, reduction = "harmony", dims = 1:30)
AIH.combined.integrated <- FindNeighbors(AIH.combined.integrated, reduction = "harmony", dims = 1:30)
AIH.combined.integrated <- FindClusters(AIH.combined.integrated, resolution = 0.4)

## had to rename these slots to get PrepSCTFindMarkers to run (see https://github.com/satijalab/seurat/issues/8235)
slot(object = AIH.combined.integrated@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Xenium"
slot(object = AIH.combined.integrated@assays$SCT@SCTModel.list[[3]], name="umi.assay")<-"Xenium"
slot(object = AIH.combined.integrated@assays$SCT@SCTModel.list[[4]], name="umi.assay")<-"Xenium"
slot(object = AIH.combined.integrated@assays$SCT@SCTModel.list[[5]], name="umi.assay")<-"Xenium"
slot(object = AIH.combined.integrated@assays$SCT@SCTModel.list[[6]], name="umi.assay")<-"Xenium"
slot(object = AIH.combined.integrated@assays$SCT@SCTModel.list[[7]], name="umi.assay")<-"Xenium"

AIH.combined.integrated <- PrepSCTFindMarkers(AIH.combined.integrated)
xenium.markers <- FindAllMarkers(AIH.combined.integrated,assay = "SCT",only.pos = TRUE)
clusters <-sort(as.numeric(as.character(unique(AIH.combined.integrated$SCT_snn_res.0.4))))
wb <- createWorkbook()
for (i in 1:length(clusters)){
  cluster <- i-1
  addWorksheet(wb,sheetName = paste0("cluster",cluster))
  writeData(wb, sheet = i, x = xenium.markers[xenium.markers$cluster==cluster,],colNames = TRUE,rowNames = FALSE)
  
}
saveWorkbook(wb,paste0(working.dir,"/results/AIH_cluster_markers_res0_4.xlsx"),overwrite = TRUE)


## Rename the clusters
# 0 - Tcells CD8
# 1 - Macrophages/Monocytes
# 2 - Fibroblasts
# 3 - Hepatocytes.central
# 4 - Kupffer cells
# 5 - Hepatocytes.mixed
# 6 - Cholangiocytes
# 7 - Tcells CD4
# 8 - Endothelial cells (vascular)
# 9 - Plasma cells
# 10 - Hepatocytes inflammatory
# 11 - Hepatocytes.portal
# 12 - dividing cells
# 13 - B-cells
# 14 - Granulocytes
# 15 - NK cells
# 16 - cDC1s
# 17 - Endothelial cells (lymphatic)
# 18 - pDCs
# 19 - Innate like T cells
# 20 - Smooth muscle cells

Idents(AIH.combined.integrated) <- "SCT_snn_res.0.4"

new.ids <- c("Tcells.CD8", #0
             
             "Mac.Monocytes", #1
             "Hepatocytes.inflammatory", #2
             "Hepatocytes.portal", #3
             "dividing.cells", #4
             "B.cells", #5
             
             "Granulocytes", #6
             "NK.cells", #7
             "cDC1s", #8
             "EC.lymphatic", #9
             "pDCs", #10
             
             "Innate.like.T.cells", #11
             "Fibroblasts", #12
             "SMC", #13
             "Hepatocytes.central", #14
             "Kupffer.cells", #15
             
             "Hepatocytes.mixed", #16
             "Cholangiocytes", #17
             "Tcells.CD4", #18
             "ECs.vascular", #19
             "Plasma.cells" #20
             )

names(new.ids) <- levels(AIH.combined.integrated)
AIH.combined.integrated <- RenameIdents(AIH.combined.integrated, new.ids)
Idents(AIH.combined.integrated) <- factor(x = Idents(AIH.combined.integrated), levels = c("Hepatocytes.central","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory","Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","SMC","Mac.Monocytes","Kupffer.cells","Granulocytes","cDC1s","pDCs","NK.cells","Innate.like.T.cells","Tcells.CD8","Tcells.CD4","B.cells","Plasma.cells", "dividing.cells"))
AIH.combined.integrated$celltype <- Idents(AIH.combined.integrated)

## 3. Add TRM1, TEM and aaCD8 annotation ---------------------------------------

AIH.combined.integrated$CD4.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="CD4",]
AIH.combined.integrated$CD69.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="CD69",]
AIH.combined.integrated$CXCR6.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="CXCR6",]
AIH.combined.integrated$IFNG.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="IFNG",]
AIH.combined.integrated$TNF.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="TNF",]
AIH.combined.integrated$SELL.expr <-  AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="SELL",]
AIH.combined.integrated$CCR7.expr <-  AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="CCR7",]

AIH.combined.integrated$TRM1.status <- AIH.combined.integrated$celltype=="Tcells.CD4"&(AIH.combined.integrated$CD69.expr>0|AIH.combined.integrated$CXCR6.expr>0)&(AIH.combined.integrated$IFNG.expr>0|AIH.combined.integrated$TNF.expr>0)
AIH.combined.integrated$aaCD8.status <- AIH.combined.integrated$celltype=="Tcells.CD8"&(AIH.combined.integrated$CD69.expr>0|AIH.combined.integrated$CXCR6.expr>0)
AIH.combined.integrated$TEM.status <- AIH.combined.integrated$celltype=="Tcells.CD4"&(AIH.combined.integrated$SELL.expr==0&AIH.combined.integrated$CD69.expr==0&AIH.combined.integrated$CCR7.expr==0)


AIH.combined.integrated$celltype.TRM1 <- as.character(AIH.combined.integrated$celltype)

AIH.combined.integrated$celltype.TRM1[AIH.combined.integrated$TRM1.status==TRUE] <- "Tcells.CD4.TRM1"
AIH.combined.integrated$celltype.TRM1[AIH.combined.integrated$aaCD8.status==TRUE] <- "Tcells.aaCD8"
AIH.combined.integrated$celltype.TRM1[AIH.combined.integrated$TEM.status==TRUE&AIH.combined.integrated$TRM1.status==FALSE] <- "Tcells.CD4.TEM"

levels.AIH.TRM1 <- c("Hepatocytes.central","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory","Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","SMC","Mac.Monocytes","Kupffer.cells","Granulocytes","cDC1s","pDCs","NK.cells","ILC","Tcells.CD8","Tcells.aaCD8","Tcells.CD4","Tcells.CD4.TRM1","Tcells.CD4.TEM","B.cells","Plasma.cells", "dividing.cells")
AIH.combined.integrated$celltype.TRM1 <- factor(AIH.combined.integrated$celltype.TRM1,levels = levels.AIH.TRM1)



saveRDS(AIH.combined.integrated,paste0(working.dir,"/SeuratObjects/AIH_integrated.rds"))
