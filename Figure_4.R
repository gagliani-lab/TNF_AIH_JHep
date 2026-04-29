## Figure 4 ##

library(Seurat) # packageVersion("Seurat") # 5.0.3
library(tidyverse) # packageVersion("tidyverse") # 2.0.0

## Insert directory here which contains the following folders:
## rawData #(contains Xenium raw data)
## results
## doc
## SeuratObjects
# working.dir <- "path/to/working/directory"

AIH.combined.integrated <- readRDS(file = paste0(working.dir,"/SeuratObjects/AIH_integrated.rds"))
Control.combined.integrated <- readRDS(file = paste0(working.dir,"/SeuratObjects/Control_integrated.rds"))

cellchat.AIH.subset <- readRDS(file = paste0(working.dir,"/SeuratObjects/AIH_CellChat.rds"))
cellchat.Control.subset <- readRDS(file = paste0(working.dir,"/SeuratObjects/Control_CellChat.rds"))


## Figure 4C -------------------------------------------------------------------

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_3__20240620__110232/")
AIH.pat2.2 <- LoadXenium(path, fov = "fov")
AIH.pat2.2 <- subset(AIH.pat2.2, subset = nCount_Xenium > 0)

AIH.pat2.2.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.2",]
AIH.pat2.2.meta <- AIH.pat2.2.meta[,colnames(AIH.pat2.2.meta) %in% c("celltype","celltype.TRM1")]
AIH.pat2.2 <- AddMetaData(AIH.pat2.2,AIH.pat2.2.meta)

AIH.pat2.2$major.celltype.3 <- as.character(AIH.pat2.2$celltype)
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype %in% c("Mac.Monocytes","Kupffer.cells","Granulocytes","cDC1s","pDCs")] <- "Myeloid.cells"
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype=="Hepatocytes.periportal.intermediate"] <- "Hepatocytes.portal"
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype=="Hepatocytes.periportal.fibrotic"] <- "Hepatocytes.mixed"
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype=="Hepatocytes.central.middle"] <- "Hepatocytes.central"
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype %in% c("Tcells.CD8","Tcells.CD4")] <- "Tcells"
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype %in% c("B.cells","Plasma.cells")] <- "Bcells"
AIH.pat2.2$major.celltype.3[AIH.pat2.2$celltype %in% c("Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","SMC", "dividing.cells","NK.cells","ILC")] <- "other"
AIH.pat2.2$major.celltype.3 <- factor(AIH.pat2.2$major.celltype.3,levels = c(
  "Hepatocytes.central","Hepatocytes.midzone.1","Hepatocytes.midzone.2","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory","Myeloid.cells","Tcells","Bcells","other"))


DefaultBoundary(AIH.pat2.2[["fov"]]) <- "segmentation"

png(paste0(working.dir,"/results/AIHpat2.2_whole_biopsy_cellborders.png"),width = 1600,height = 3200,res = 300)
ImageDimPlot(AIH.pat2.2, fov = "fov", axes = TRUE, border.color = "darkgrey", border.size = 0, nmols = 10000,group.by = "major.celltype.3",dark.background = FALSE,coord.fixed = TRUE, crop = TRUE) +
  scale_fill_manual(values=c(rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5")),"#A1D99B","#E31A1C","#FED976","#969696"),drop = FALSE) +
  theme_classic()
dev.off()

png(paste0(working.dir,"/results/AIHpat2.2_whole_biopsy_cellborders_detail.png"),width = 1600,height = 1600,res = 300)
ImageDimPlot(AIH.pat2.2, fov = "fov", axes = TRUE, border.color = "darkgrey", border.size = 0, nmols = 10000,group.by = "major.celltype.3",dark.background = FALSE,coord.fixed = TRUE, crop = TRUE) +
  scale_fill_manual(values=c(rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5")),"#A1D99B","#E31A1C","#FED976","#969696"),drop = FALSE) +
  theme_classic() + 
  scale_y_continuous(limits = c(2000,4000)) + 
  scale_x_continuous(limits = c(1000,2000))
dev.off()


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

DefaultBoundary(AIH.pat3[["fov"]]) <- "segmentation"
png(paste0(working.dir,"/results/AIHpat3_whole_biopsy_cellborders.png"),width = 1600,height = 3200,res = 300)
ImageDimPlot(AIH.pat3, fov = "fov", axes = TRUE, border.color = "darkgrey", border.size = 0, nmols = 10000,group.by = "major.celltype.3",dark.background = FALSE,coord.fixed = TRUE, crop = TRUE) +
  scale_fill_manual(values=c(rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5")),"#A1D99B","#E31A1C","#FED976","#969696"),drop = FALSE) +
  theme_classic()
dev.off()

png(paste0(working.dir,"/results/AIHpat3_whole_biopsy_cellborders_detail.png"),width = 1600,height = 1600,res = 300)
ImageDimPlot(AIH.pat3, fov = "fov", axes = TRUE, border.color = "darkgrey", border.size = 0, nmols = 10000,group.by = "major.celltype.3",dark.background = FALSE,coord.fixed = TRUE, crop = TRUE) +
  scale_fill_manual(values=c(rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5")),"#A1D99B","#E31A1C","#FED976","#969696"),drop = FALSE) +
  theme_classic() + 
  scale_y_continuous(limits = c(7800,9400)) + 
  scale_x_continuous(limits = c(1200,2900))
dev.off()


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
Control.pat3.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat3",]
Control.pat3.meta <- Control.pat3.meta[,colnames(Control.pat3.meta) %in% c("celltype","celltype.TRM1")]
Control.pat3 <- AddMetaData(Control.pat3,Control.pat3.meta)

Control.pat3$major.celltype.3 <- as.character(Control.pat3$celltype)
Control.pat3$major.celltype.3[Control.pat3$celltype %in% c("Mac.Monocytes","Kupffer.cells","Granulocytes","cDC1s","pDCs")] <- "Myeloid.cells"
Control.pat3$major.celltype.3[Control.pat3$celltype=="Hepatocytes.periportal"] <- "Hepatocytes.portal"
Control.pat3$major.celltype.3[Control.pat3$celltype=="Hepatocytes.central.middle"] <- "Hepatocytes.midzone.1"
Control.pat3$major.celltype.3[Control.pat3$celltype=="Hepatocytes.periportal.middle"] <- "Hepatocytes.midzone.2"
Control.pat3$major.celltype.3[Control.pat3$celltype %in% c("Tcells.CD8","Tcells.CD4")] <- "Tcells"
Control.pat3$major.celltype.3[Control.pat3$celltype %in% c("B.cells","Plasma.cells")] <- "Bcells"
Control.pat3$major.celltype.3[Control.pat3$celltype %in% c("Cholangiocytes","Fibroblasts","ECs.vascular","EC.lymphatic","SMC", "dividing.cells","NK.cells","ILC","Neutrophils","undefined")] <- "other"
Control.pat3$major.celltype.3 <- factor(Control.pat3$major.celltype.3,levels = c(
  "Hepatocytes.central","Hepatocytes.midzone.1","Hepatocytes.midzone.2","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory","Myeloid.cells","Tcells","Bcells","other"))

DefaultBoundary(Control.pat3[["fov"]]) <- "segmentation"

png(paste0(working.dir,"/results/Controlpat3_whole_biopsy_cellborders.png"),width = 1600,height = 3200,res = 300)
ImageDimPlot(Control.pat3, fov = "fov", axes = TRUE, border.color = "darkgrey", border.size = 0, nmols = 10000,group.by = "major.celltype.3",dark.background = FALSE,coord.fixed = TRUE, crop = TRUE) +
  scale_fill_manual(values=c(rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5")),"#A1D99B","#E31A1C","#FED976","#969696"),drop = FALSE) +
  theme_classic()
dev.off()

png(paste0(working.dir,"/results/Controlpat3_whole_biopsy_cellborders_detail.png"),width = 1600,height = 1600,res = 300)
ImageDimPlot(Control.pat3, fov = "fov", axes = TRUE, border.color = "darkgrey", border.size = 0, nmols = 10000,group.by = "major.celltype.3",dark.background = FALSE,coord.fixed = TRUE, crop = TRUE) +
  scale_fill_manual(values=c(rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5")),"#A1D99B","#E31A1C","#FED976","#969696"),drop = FALSE) +
  theme_classic() + 
  scale_y_continuous(limits = c(1300,2500)) + 
  scale_x_continuous(limits = c(5200,6100))
dev.off()

## Figure 4D -------------------------------------------------------------------

cropped.coords.A <- Crop(AIH.pat3[["fov"]], y = c(8150, 8450), x = c(1600, 1800), coords = "plot")
cropped.coords.B <- Crop(AIH.pat3[["fov"]], y = c(1500, 1750), x = c(800, 1200), coords = "plot")

AIH.pat3[["sectionA"]] <- cropped.coords.A
AIH.pat3[["sectionB"]] <- cropped.coords.B

DefaultBoundary(AIH.pat3[["sectionA"]]) <- "segmentation"
DefaultBoundary(AIH.pat3[["sectionB"]]) <- "segmentation"

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

pdf(paste0(working.dir,"/results/AIHpat3_detailA_no_dots.pdf"),width = 15,height = 9)
ImageDimPlot(AIH.pat3, fov = "sectionA", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red, nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,coord.fixed = TRUE) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailA_IL15.pdf"),width = 15,height = 9)
ImageDimPlot(AIH.pat3, fov = "sectionA", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("IL15","IL15RA","IL2RB"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1","turquoise3"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailA_TNF.pdf"),width = 15,height = 9)
ImageDimPlot(AIH.pat3, fov = "sectionA", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("TNF","TNFRSF1A","TNFRSF1B"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1","turquoise3"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailA_Granzymes.pdf"),width = 15,height = 9)
ImageDimPlot(AIH.pat3, fov = "sectionA", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("GZMA","GZMB","GZMK","GZMH","GZMM"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1","turquoise3","gold","darkgreen"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailA_IL7_IL7R.pdf"),width = 15,height = 9)
ImageDimPlot(AIH.pat3, fov = "sectionA", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("IL7","IL7R"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1"),mols.size = 0.5,coord.fixed = TRUE) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

### Detail B

my.ratio=(375/250)

pdf(paste0(working.dir,"/results/AIHpat3_detailB_no_dots.pdf"),width = 15,height = 15)
ImageDimPlot(AIH.pat3, fov = "sectionB", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red, nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,coord.fixed = TRUE) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailB_IL15.pdf"),width = 15,height = 15)
ImageDimPlot(AIH.pat3, fov = "sectionB", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("IL15","IL15RA","IL2RB"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1","turquoise3"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailB_TNF.pdf"),width = 15,height = 15)
ImageDimPlot(AIH.pat3, fov = "sectionB", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("TNF","TNFRSF1A","TNFRSF1B"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1","turquoise3"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailB_Granzymes.pdf"),width = 15,height = 15)
ImageDimPlot(AIH.pat3, fov = "sectionB", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("GZMA","GZMB","GZMK","GZMH","GZMM"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = FALSE,mols.cols = c("deeppink","purple1","turquoise3","gold","darkgreen"),mols.size = 0.5) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()

pdf(paste0(working.dir,"/results/AIHpat3_detailB_IL7_IL7R.pdf"),width = 15,height = 15)
ImageDimPlot(AIH.pat3, fov = "sectionB", axes = TRUE, border.color = "darkgrey", border.size = 0.2, cols = colors.AIH.TRM1.red,molecules=c("IL7","IL7R"), nmols = 10000,group.by = "celltype.TRM1.simple",dark.background = TRUE,mols.cols = c("deeppink","purple1"),mols.size = 0.5,coord.fixed = TRUE) & scale_fill_manual(values = colors.AIH.TRM1.red,drop=FALSE) & theme_classic() & coord_flip() & theme(aspect.ratio = my.ratio)
dev.off()


## Figure 4E -------------------------------------------------------------------

cytos <- c("TNF", "IL15", "IL12AB", "IL7", "IFNG", "IL18",  "FASLG","SELE","CDH1","ICAM1","IL1B")
CellChatDB.use <- subsetDB(CellChatDB.human, 
                           search = cytos, 
                           key = c("ligand"))

# Add IL15+IL15RA -> IL2RG/Il2RB 
df.add <- data.frame(colnames(CellChatDB.use$interaction)) %>% t %>% as.data.frame
df.add[1,] <- rep("", times = ncol(df.add))
colnames(df.add) <- colnames(CellChatDB.use$interaction)
df.add$pathway_name = "IL2"
df.add$interaction_name = "IL15RA_IL2RG_IL2RB"
df.add$interaction_name_2 = "IL15RA - (IL2RG+IL2RB)"
df.add$ligand = "IL15RA"
df.add$ligand.symbol = "IL15RA"
df.add$receptor = "IL2RG_IL2RB"
df.add$receptor.symbol = "IL2RG, IL2RB"
rownames(df.add) = "IL15RA_IL2RG_IL2RB"


# Add IL15RA + IL15 -> IL2RG/Il2RB 
df.add2 <- data.frame(colnames(CellChatDB.use$interaction)) %>% t %>% as.data.frame
df.add2[1,] <- rep("", times = ncol(df.add2))
colnames(df.add2) <- colnames(CellChatDB.use$interaction)
df.add2$pathway_name = "IL2"
df.add2$interaction_name = "IL15RA_IL15_IL2RG_IL2RB"
df.add2$interaction_name_2 = "(IL15RA+IL15) - (IL2RG+IL2RB)"
df.add2$ligand = "IL15RA_IL15"
df.add2$ligand.symbol = "IL15RA, IL15"
df.add2$receptor = "IL2RG_IL2RB"
df.add2$receptor.symbol = "IL2RG, IL2RB"
rownames(df.add2) = "IL15RA_IL15_IL2RG_IL2RB"

# Add Complex Entry
df.complex <- data.frame(subunit_1 = c("IL2RG", "IL15RA"), 
                         subunit_2 = c("IL2RB", "IL15"), 
                         subunit_3 = c("", ""), 
                         subunit_4 =  c("", ""), subunit_5 =  c("", ""), row.names = c("IL2RG_IL2RB", "IL15RA_IL15"))

CellChatDB.use$interaction <- rbind(CellChatDB.use$interaction, df.add,df.add2)
CellChatDB.use$complex <- rbind(CellChatDB.use$complex, df.complex)


cellchat.AIH.subset@DB <- CellChatDB.use
cellchat.AIH.subset <- subsetData(cellchat.AIH.subset)

# We want to include all selected ligand -> no DE, no threshold for filtering
cellchat.AIH.subset <- identifyOverExpressedGenes(cellchat.AIH.subset,do.fast = T, 
                                                     do.DE = FALSE, thresh.p = 1, thresh.pc = 0, 
                                                     thresh.fc = 0)
cellchat.AIH.subset <- identifyOverExpressedInteractions(cellchat.AIH.subset)

# Threshold trim 0.01, mainly because of Il15RA, which was expressed by only 4.7% of aa CD8 cells 
cellchat.AIH.subset <- computeCommunProb(cellchat.AIH.subset, type="truncatedMean", trim = 0.01,distance.use = FALSE, interaction.range = 50, scale.distance = NULL, contact.dependent = TRUE, contact.range = 10, population.size = TRUE)

cellchat.AIH.subset <- filterCommunication(cellchat.AIH.subset, min.cells = 10)

cellchat.AIH.subset <- computeCommunProbPathway(cellchat.AIH.subset, thresh = 1)
cellchat.AIH.subset <- aggregateNet(cellchat.AIH.subset)


cellchat.Control.subset@DB <- CellChatDB.use
cellchat.Control.subset <- subsetData(cellchat.Control.subset)


# We want to include all selected ligand -> no DE, no threshold for filtering
cellchat.Control.subset <- identifyOverExpressedGenes(cellchat.Control.subset,do.fast = T, 
                                                         do.DE = FALSE, thresh.p = 1, thresh.pc = 0, 
                                                         thresh.fc = 0)
cellchat.Control.subset <- identifyOverExpressedInteractions(cellchat.Control.subset)

# Threshold trim 0.01, mainly because of Il15RA, which was expressed by only 4.7% of aa CD8 cells 
cellchat.Control.subset <- computeCommunProb(cellchat.Control.subset, type="truncatedMean", trim = 0.01,distance.use = FALSE, interaction.range = 50, scale.distance = NULL, contact.dependent = TRUE, contact.range = 10, population.size = TRUE)

cellchat.Control.subset <- filterCommunication(cellchat.Control.subset, min.cells = 10)

cellchat.Control.subset <- computeCommunProbPathway(cellchat.Control.subset, thresh = 1)
cellchat.Control.subset <- aggregateNet(cellchat.Control.subset)

object.list = list(cellchat.Control.subset,cellchat.AIH.subset)
names(object.list) = c("Control","AIH")
cellchat.all <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# get only significant interactions 
df.net <- subsetCommunication(cellchat.all, thresh = 0.05)
df.net$Control$disease <- "Control"
df.net$AIH$disease <- "AIH"
df.net <- rbind(df.net$Control,df.net$AIH)

df.net$source <- factor(df.net$source,levels = c("Hepatocytes","Myeloid.cells","Tcells.CD4.TRM1","Tcells.aaCD8"))
df.net$target <- factor(df.net$target,levels = c("Hepatocytes","Myeloid.cells","Tcells.CD4.TRM1","Tcells.aaCD8"))

df.net.AIH <- df.net[df.net$disease=="AIH",]
df.net.Control <- df.net[df.net$disease=="Control",]

df.net.AIH.red <- df.net.AIH[df.net.AIH$interaction_name %in% c("IFNG_IFNGR1_IFNGR2","TNF_TNFRSF1A","TNF_TNFRSF1B","FASL_FAS","IL7_IL7R_IL2RG","IL15_IL15RA_IL2RB","IL18 - (IL18R1+IL18RAP)","IL15RA_IL2RG_IL2RB","IL15RA_IL15_IL2RG_IL2RB"),]
df.net.AIH.red$interaction_name_2 <- factor(df.net.AIH.red$interaction_name_2,levels = rev(c("TNF - TNFRSF1B","TNF - TNFRSF1A","IFNG - (IFNGR1+IFNGR2)","FASL - FAS","IL7 - (IL7R+IL2RG)","IL18 - (IL18R1+IL18RAP)","IL15 - (IL15RA+IL2RB)","IL15RA - (IL2RG+IL2RB)","(IL15RA+IL15) - (IL2RG+IL2RB)")))


df.net.AIH.red <- df.net.AIH.red[!((df.net.AIH.red$source=="Myeloid.cells"&df.net.AIH.red$target %in% c("Hepatocytes","Myeloid.cells"))|(df.net.AIH.red$source=="Hepatocytes"&df.net.AIH.red$target %in% c("Hepatocytes","Myeloid.cells"))),]

df.net.Control.red <- df.net.Control[df.net.Control$interaction_name %in% c("IFNG_IFNGR1_IFNGR2","TNF_TNFRSF1A","TNF_TNFRSF1B","FASL_FAS","IL7_IL7R_IL2RG","IL15_IL15RA_IL2RB","IL15RA_IL2RG_IL2RB","IL15RA_IL15_IL2RG_IL2RB"),]
df.net.Control.red$interaction_name_2 <- factor(df.net.Control.red$interaction_name_2,levels = rev(c("TNF - TNFRSF1B","TNF - TNFRSF1A","IFNG - (IFNGR1+IFNGR2)","FASL - FAS","IL7 - (IL7R+IL2RG)","IL18 - (IL18R1+IL18RAP)","IL15 - (IL15RA+IL2RB)","IL15RA - (IL2RG+IL2RB)","(IL15RA+IL15) - (IL2RG+IL2RB)")))


df.net.Control.red <- df.net.Control.red[!((df.net.Control.red$source=="Myeloid.cells"&df.net.Control.red$target %in% c("Hepatocytes","Myeloid.cells"))|(df.net.Control.red$source=="Hepatocytes"&df.net.Control.red$target %in% c("Hepatocytes","Myeloid.cells"))),]

interactions <- c("IFNG_IFNGR1_IFNGR2","TNF_TNFRSF1A","TNF_TNFRSF1B","FASL_FAS","IL7_IL7R_IL2RG","IL15_IL15RA_IL2RB","IL15RA_IL2RG_IL2RB","IL15RA_IL15_IL2RG_IL2RB")

pdf(paste0(working.dir,"/results/CC_all_heatmaps_grouped.pdf"),width = 9,height = 4.5)

for (i in 1:length(interactions)){
  df.net.Control.red <- df.net.Control[df.net.Control$interaction_name==interactions[i],]
  df.net.Control.red.summed <- df.net.Control.red %>%
    group_by(interaction_name) %>%
    mutate(st = paste0(source, " > ", target))
  
  df.net.AIH.red <- df.net.AIH[df.net.AIH$interaction_name==interactions[i],]
  df.net.AIH.red.summed <- df.net.AIH.red %>%
    group_by(interaction_name) %>%
    mutate(st = paste0(source, " > ", target))
  df.net.all.red.summed <- rbind(df.net.AIH.red.summed,df.net.Control.red.summed)
  df.net.all.red.summed$disase <- factor(df.net.all.red.summed$disease,levels = c("AIH","Control"))
  
  p <- ggplot(df.net.all.red.summed,aes(y=disease, x=interaction(target, source, sep="&"), fill=prob)) +
    geom_tile() +
    scale_x_discrete(guide = guide_axis_nested(delim = "&"),drop = FALSE) + 
    scale_y_discrete(drop=FALSE) +
    colors$palette$DarkBlue_Red$ggfill() + 
    labs(x="", y="") + 
    theme_paper(aspect.ratio = 1/8) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    labs(title=df.net.all.red.summed$interaction_name_2[1]) +
    coord_equal()
  
  options(repr.plot.height=7, repr.plot.width=7)
  
  print(p)
}
dev.off()