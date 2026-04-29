## Extended Figure 6 ##

library(Seurat) # packageVersion("Seurat") # 5.0.3
library(tidyverse) # packageVersion("tidyverse") # 2.0.0

## Insert directory here which contains the following folders:
## rawData #(contains Xenium raw data)
## results
## doc
## SeuratObjects
# working.dir <- "path/to/working/directory"

AIH.combined.integrated <- readRDS(file = paste0(working.dir,"/SeuratObjects/AIH_integrated.rds"))

## Extended Figure 8A - Xenium detail -----------------------------------------------

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_4__20240620__110232/")
AIH.pat3 <- LoadXenium(path, fov = "fov")
AIH.pat3 <- subset(AIH.pat3, subset = nCount_Xenium > 0)
AIH.pat3.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat3",]
AIH.pat3.meta <- AIH.pat3.meta[,colnames(AIH.pat3.meta) %in% c("celltype","celltype.TRM1")]
AIH.pat3 <- AddMetaData(AIH.pat3,AIH.pat3.meta)

AIH.pat3$major.celltype.3 <- as.character(AIH.pat3$celltype)
AIH.pat3$major.celltype.3[AIH.pat3$celltype %in% c("Mac.Monocytes","Kupffer.cells","Granulocytes","cDC1s","pDCs")] <- "Myeloid.cells"
AIH.pat3$major.celltype.3[AIH.pat3$celltype=="Hepatocytes.periportal.intermediate"] <- "Hepatocytes.portal"
AIH.pat3$major.celltype.3[AIH.pat3$celltype=="Hepatocytes.periportal.fibrotic"] <- "Hepatocytes.mixed"
AIH.pat3$major.celltype.3[AIH.pat3$celltype=="Hepatocytes.central.middle"] <- "Hepatocytes.central"
AIH.pat3$major.celltype.3[AIH.pat3$celltype %in% c("Tcells.CD8","Tcells.CD4")] <- "Tcells"
AIH.pat3$major.celltype.3[AIH.pat3$celltype %in% c("B.cells","Plasma.cells")] <- "Bcells"
AIH.pat3$major.celltype.3[AIH.pat3$celltype %in% c("Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","SMC", "dividing.cells","NK.cells","ILC")] <- "other"
AIH.pat3$major.celltype.3 <- factor(AIH.pat3$major.celltype.3,levels = c(
  "Hepatocytes.central","Hepatocytes.midzone.1","Hepatocytes.midzone.2","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory","Myeloid.cells","Tcells","Bcells","other"))

cropped.coords.A <- Crop(AIH.pat3[["fov"]], y = c(8150, 8450), x = c(1600, 1800), coords = "plot")

AIH.pat3[["sectionA"]] <- cropped.coords.A

DefaultBoundary(AIH.pat3[["sectionA"]]) <- "segmentation"

AIH.pat3$celltype.TRM1.simple <- rep("other",ncol(AIH.pat3))

AIH.pat3$celltype.TRM1.simple[grepl("Hepat",AIH.pat3$celltype.TRM1)] <- "Hepatocytes"
AIH.pat3$celltype.TRM1.simple[AIH.pat3$celltype.TRM1=="Fibroblasts"] <- "Fibroblasts"

AIH.pat3$celltype.TRM1.simple[AIH.pat3$celltype.TRM1=="Cholangiocytes"] <- "Cholangiocytes"
AIH.pat3$celltype.TRM1.simple[grepl("EC",AIH.pat3$celltype.TRM1)] <- "Endothelial.cells"
AIH.pat3$celltype.TRM1.simple[AIH.pat3$celltype.TRM1=="Tcells.CD8"] <- "Tcells.CD8"
AIH.pat3$celltype.TRM1.simple[AIH.pat3$celltype.TRM1=="Tcells.aaCD8"] <- "Tcells.aaCD8"
AIH.pat3$celltype.TRM1.simple[AIH.pat3$celltype.TRM1=="Tcells.CD4.TRM1"] <- "Tcells.CD4.TRM1"

AIH.pat3$celltype.TRM1.simple[AIH.pat3$celltype.TRM1 %in% c("Mac.Monocytes","Kupffer.cells","Granulocytes","pDCs","cDC1s")] <- "Myeloid.cells"

AIH.pat3$celltype.TRM1.simple <- factor(AIH.pat3$celltype.TRM1.simple,levels = c("Hepatocytes","Cholangiocytes","Fibroblasts","Endothelial.cells","Tcells.CD8","Tcells.aaCD8","Tcells.CD4.TRM1","Myeloid.cells","other"))

## Detail A

my.ratio=(200/275)

pdf(paste0(working.dir,"/results/AIHpat3_detailA_ITGAL_ICAM1.pdf"),width = 15,height = 9)
ImageDimPlot(AIH.pat3, fov = "sectionA", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("ITGAL","ICAM1"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()