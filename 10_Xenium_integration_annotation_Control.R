## Integration + Annotation Xenium Control data ##

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

## 1. Import Control data ------------------------------------------------------

### 1.1 Control Pat 2

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_6__20240620__110232/")
Control.pat2 <- LoadXenium(path, fov = "fov")
Control.pat2 <- subset(Control.pat2, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(Control.pat2),]
csv.cells <- csv.cells[match(colnames(Control.pat2),csv.cells$cell_id),]
all(colnames(Control.pat2)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

Control.pat2 <- AddMetaData(Control.pat2,csv.cells)
Control.pat2 <- RenameCells(Control.pat2,add.cell.id = "Control.pat2")
Control.pat2$sample <- "Control.pat2"


### 1.2 Control Pat 3
# Cut this so only the piece at the bottom is included (the top pieces are from an AIH biopsy)

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_3__20240620__110232/")
Control.pat3 <- LoadXenium(path, fov = "fov")
Control.pat3 <- subset(Control.pat3, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(Control.pat3),]
csv.cells <- csv.cells[match(colnames(Control.pat3),csv.cells$cell_id),]
all(colnames(Control.pat3)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

Control.pat3 <- AddMetaData(Control.pat3,csv.cells)
Control.pat3 <- subset(Control.pat3, subset = y_centroid>4500)
Control.pat3 <- RenameCells(Control.pat3,add.cell.id = "Control.pat3")
Control.pat3$sample <- "Control.pat3"


### 1.3 Control Pat 5

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_5__20240620__110232/")
Control.pat5 <- LoadXenium(path, fov = "fov")
Control.pat5 <- subset(Control.pat5, subset = nCount_Xenium > 0)
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(Control.pat5),]
csv.cells <- csv.cells[match(colnames(Control.pat5),csv.cells$cell_id),]
all(colnames(Control.pat5)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

Control.pat5 <- AddMetaData(Control.pat5,csv.cells)
Control.pat5 <- RenameCells(Control.pat5,add.cell.id = "Control.pat5")
Control.pat5$sample <- "Control.pat5"

## 2. Integrate Control data ---------------------------------------------------

Control.combined.integrated <- merge(Control.pat2,y = c(Control.pat3,Control.pat5))
Control.combined.integrated <- SCTransform(Control.combined.integrated, assay = "Xenium")
Control.combined.integrated <- RunPCA(Control.combined.integrated, npcs = 30, features = rownames(Control.combined.integrated))
Control.combined.integrated <- Control.combined.integrated %>% RunHarmony("sample",plot_convergence = TRUE,max_iter = 15,early_stop=FALSE)
Control.combined.integrated <- RunUMAP(Control.combined.integrated, reduction = "harmony", dims = 1:30)
Control.combined.integrated <- FindNeighbors(Control.combined.integrated, reduction = "harmony", dims = 1:30)
Control.combined.integrated <- FindClusters(Control.combined.integrated, resolution = 0.4)

### 2.1 Subcluster the T-cells

Idents(Control.combined.integrated) <- "SCT_snn_res.0.4"
Control.combined.integrated <- FindSubCluster(Control.combined.integrated,7,graph.name = "SCT_snn",resolution = 0.3)
Idents(Control.combined.integrated) <- "sub.cluster"
Control.combined.integrated$sub.cluster <- factor(Control.combined.integrated$sub.cluster, levels = c("0","1","2","3","4","5","6","7_0","7_1","7_2","7_3","7_4","7_5","7_6","8","9","10","11","12","13","14"))

Idents(Control.combined.integrated) <- "sub.cluster"
## had to rename these slots to get PrepSCTFindMarkers to run (see https://github.com/satijalab/seurat/issues/8235)
slot(object = Control.combined.integrated@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Xenium"
slot(object = Control.combined.integrated@assays$SCT@SCTModel.list[[3]], name="umi.assay")<-"Xenium"

Control.combined.integrated <- PrepSCTFindMarkers(Control.combined.integrated)
xenium.markers <- FindAllMarkers(Control.combined.integrated,assay = "SCT",only.pos = TRUE)
clusters <- levels(Control.combined.integrated$sub.cluster)

wb <- createWorkbook()
for (i in 1:length(clusters)){
  cluster <- clusters[i]
  addWorksheet(wb,sheetName = paste0("cluster",cluster))
  writeData(wb, sheet = i, x = xenium.markers[xenium.markers$cluster==cluster,],colNames = TRUE,rowNames = FALSE)
  
}
saveWorkbook(wb,paste0(working.dir,"/results/Control_cluster_markers_res0_4.xlsx"),overwrite = TRUE)

## Rename the clusters
# 0 - Hepatocytes.midzone.1
# 1 - Hepatocytes.midzone.2
# 2 - Hepatocytes.portal
# 3 - Endothelial cells (vascular)
# 4 - Fibroblasts
# 5 - Kupffer cells
# 6 - Hepatocytes central
# 7_0 - Tcells CD8
# 7_1 - NK cells
# 7_2 - NK cells
# 7_3 - undefined
# 7_4 - Tcells CD4
# 7_5 - undefined
# 7_6 - undefined
# 8 - Macrophages/Monocytes
# 9 - B-cells
# 10 - Endothelial cells (lymphatic)
# 11 - Cholangiocytes
# 12 - Fibroblasts
# 13 - Plasma cells
# 14 - Hepatocytes inflammatory


Idents(Control.combined.integrated) <- "sub.cluster"

new.ids <- c("Hepatocytes.midzone.1", #0
             
             "Hepatocytes.midzone.2", #1
             "Hepatocytes.portal", #2
             "ECs.vascular", #3
             "Fibroblasts", #4
             "Kupffer.cells", #5
             
             "Hepatocytes.central", #6
             "Tcells.CD8", #7_0
             "NK.cells", #7_1
             "NK.cells", #7_2
             "undefined", #7_3
             
             "Tcells.CD4", #7_4
             "undefined", #7_5
             "undefined", #7_6
             "Mac.Monocytes", #8
             "B.cells", #9
             "EC.lymphatic", #10
             
             "Cholangiocytes", #11
             "Fibroblasts", #12
             "Plasma.cells", #13
             "Hepatocytes.inflammatory" #14
             )

names(new.ids) <- levels(Control.combined.integrated)
Control.combined.integrated <- RenameIdents(Control.combined.integrated, new.ids)
Idents(Control.combined.integrated) <- factor(x = Idents(Control.combined.integrated), levels = c("Hepatocytes.central","Hepatocytes.midzone.1","Hepatocytes.midzone.2","Hepatocytes.portal","Hepatocytes.inflammatory","Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","Mac.Monocytes","Kupffer.cells","NK.cells","Tcells.CD8","Tcells.CD4","B.cells","Plasma.cells", "undefined"))
Control.combined.integrated$celltype <- Idents(Control.combined.integrated)

## 3. Add TRM1, TEM and aaCD8 annotation ---------------------------------------

Control.combined.integrated$CD4.expr <- Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="CD4",]
Control.combined.integrated$CD69.expr <- Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="CD69",]
Control.combined.integrated$CXCR6.expr <- Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="CXCR6",]
Control.combined.integrated$IFNG.expr <- Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="IFNG",]
Control.combined.integrated$TNF.expr <- Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="TNF",]
Control.combined.integrated$SELL.expr <-  Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="SELL",]
Control.combined.integrated$CCR7.expr <-  Control.combined.integrated[["Xenium"]]$counts[rownames(Control.combined.integrated[["Xenium"]]$counts)=="CCR7",]

Control.combined.integrated$TRM1.status <- Control.combined.integrated$celltype=="Tcells.CD4"&(Control.combined.integrated$CD69.expr>0|Control.combined.integrated$CXCR6.expr>0)&(Control.combined.integrated$IFNG.expr>0|Control.combined.integrated$TNF.expr>0)

Control.combined.integrated$aaCD8.status <- Control.combined.integrated$celltype=="Tcells.CD8"&(Control.combined.integrated$CD69.expr>0|Control.combined.integrated$CXCR6.expr>0)

Control.combined.integrated$TEM.status <- Control.combined.integrated$celltype=="Tcells.CD4"&(Control.combined.integrated$SELL.expr==0&Control.combined.integrated$CD69.expr==0&Control.combined.integrated$CCR7.expr==0)


Control.combined.integrated$celltype.TRM1 <- as.character(Control.combined.integrated$celltype)

Control.combined.integrated$celltype.TRM1[Control.combined.integrated$TRM1.status==TRUE] <- "Tcells.CD4.TRM1"
Control.combined.integrated$celltype.TRM1[Control.combined.integrated$aaCD8.status==TRUE] <- "Tcells.aaCD8"
Control.combined.integrated$celltype.TRM1[Control.combined.integrated$TEM.status==TRUE&Control.combined.integrated$TRM1.status==FALSE] <- "Tcells.CD4.TEM"

levels.Control.TRM1 <- c("Hepatocytes.central","Hepatocytes.midzone.1","Hepatocytes.midzone.2","Hepatocytes.portal","Hepatocytes.inflammatory","Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","Mac.Monocytes","Kupffer.cells","NK.cells","Tcells.CD8","Tcells.aaCD8","Tcells.CD4","Tcells.CD4.TRM1","Tcells.CD4.TEM","B.cells","Plasma.cells", "undefined")
Control.combined.integrated$celltype.TRM1 <- factor(Control.combined.integrated$celltype.TRM1,levels = levels.Control.TRM1)


saveRDS(Control.combined.integrated,paste0(working.dir,"/SeuratObjects/Control_integrated.rds"))


