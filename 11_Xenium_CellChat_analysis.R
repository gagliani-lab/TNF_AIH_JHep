## Integration + Annotation Xenium Control data ##

library(Seurat) # packageVersion("Seurat") # 5.0.3
library(tidyverse) # packageVersion("tidyverse") # 2.0.0
library(ggh4x) # packageVersion("ggh4x") # 0.2.8
library(CellChat) # packageVersion("CellChat") # 2.1.2


## Insert directory here which contains the following folders:
## rawData #(contains Xenium raw data)
## results
## doc
## SeuratObjects
# working.dir <- "path/to/working/directory"

colors <- list()
colors$palette$DarkBlue_Red$ggfill <- function(midpoint = 0)  scale_fill_gradient2(low = "#3F517C",high = "#EB5D12",mid = "whitesmoke",midpoint = midpoint,na.value = "white")

theme_paper <- 
  function(bordertype = "closed", # closed (surrounded by border), open (x an y axis), blank (no axis line / border)
           aspect.ratio = 1/1, # ratio of plot -> y/x
           base_size = 12, # font size
           base_family = "Helvetica", # font familiy
           legend_position = "right", 
           title_size= base_size){
    
    theme_closed <- 
      theme_bw(base_size = base_size, base_family = base_family) %+replace%
      
      theme(
        aspect.ratio = aspect.ratio, 
        # Panel
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "#000000",
                                    linewidth= unit(units = "npc", x = 0.5)), 
        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # Axis
        axis.line = element_blank(), axis.ticks = element_line(colour = "#000000"),
        axis.text = element_text(size = base_size), plot.title=element_text(size = title_size),
        axis.title = element_text(size = base_size), legend.text=element_text(size = base_size),
        # Legend
        legend.key = element_blank(),
        legend.title=element_text(size = base_size), strip.text=element_text(size = base_size),
        legend.position = legend_position, 
        legend.background = element_blank()
      )
    
    if(bordertype == "open")
      return(theme_closed %+replace% theme(axis.line = element_line(colour = "#000000", linewidth= unit(units = "npc", x = 0.5)), panel.border = element_blank()))
    else if(bordertype == "blank")
      return(theme_closed %+replace% theme(axis.line = element_blank(), panel.border = element_blank()))
    
    return(theme_closed)
    
  }


AIH.combined.integrated <- readRDS(paste0(working.dir,"/SeuratObjects/AIH_integrated.rds"))
Control.combined.integrated <- readRDS(paste0(working.dir,"/SeuratObjects/Control_integrated.rds"))


## 1. Load objects and annotation

### 1.1 Load AIH biopsies and filter out anything labeled "other"
## Keep only Hepatocytes, aaCD8, TRM1, Myeloid cells

#### 1.1.1 AIH Pat 1
path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_1__20240620__110232/")
AIH.pat1 <- LoadXenium(path, fov = "fov")
AIH.pat1 <- subset(AIH.pat1, subset = nCount_Xenium > 0)
AIH.pat1 <- SCTransform(AIH.pat1, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat1),]
csv.cells <- csv.cells[match(colnames(AIH.pat1),csv.cells$cell_id),]
all(colnames(AIH.pat1)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat1 <- AddMetaData(AIH.pat1,csv.cells)
AIH.pat1.meta <- as.data.frame(AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat1",])
AIH.pat1.meta <- AIH.pat1.meta[,colnames(AIH.pat1.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat1 <- AddMetaData(AIH.pat1,AIH.pat1.meta)
AIH.pat1 <- subset(AIH.pat1,subset = celltype.simple!="other")
AIH.pat1 <- RenameCells(AIH.pat1,add.cell.id = "AIH.pat1")
AIH.pat1$sample <- "AIH.pat1"


#### 1.1.2 AIH Pat 2 Block 1

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_2__20240620__110232/")
AIH.pat2.1 <- LoadXenium(path, fov = "fov")
AIH.pat2.1 <- subset(AIH.pat2.1, subset = nCount_Xenium > 0)
AIH.pat2.1 <- SCTransform(AIH.pat2.1, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat2.1),]
csv.cells <- csv.cells[match(colnames(AIH.pat2.1),csv.cells$cell_id),]
all(colnames(AIH.pat2.1)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat2.1 <- AddMetaData(AIH.pat2.1,csv.cells)
AIH.pat2.1.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.1",]
AIH.pat2.1.meta <- AIH.pat2.1.meta[,colnames(AIH.pat2.1.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat2.1 <- AddMetaData(AIH.pat2.1,AIH.pat2.1.meta)
AIH.pat2.1 <- subset(AIH.pat2.1,subset = celltype.simple!="other")
AIH.pat2.1 <- RenameCells(AIH.pat2.1,add.cell.id = "AIH.pat2.1")
AIH.pat2.1$sample <- "AIH.pat2.1"


#### 1.1.3 AIH Pat 2 Block 2

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_3__20240620__110232/")
AIH.pat2.2 <- LoadXenium(path, fov = "fov")
AIH.pat2.2 <- subset(AIH.pat2.2, subset = nCount_Xenium > 0)
AIH.pat2.2 <- SCTransform(AIH.pat2.2, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat2.2),]
csv.cells <- csv.cells[match(colnames(AIH.pat2.2),csv.cells$cell_id),]
all(colnames(AIH.pat2.2)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat2.2 <- AddMetaData(AIH.pat2.2,csv.cells)
AIH.pat2.2.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.2",]
AIH.pat2.2.meta <- AIH.pat2.2.meta[,colnames(AIH.pat2.2.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat2.2 <- AddMetaData(AIH.pat2.2,AIH.pat2.2.meta)
AIH.pat2.2 <- subset(AIH.pat2.2,subset = celltype.simple!="other")
AIH.pat2.2 <- RenameCells(AIH.pat2.2,add.cell.id = "AIH.pat2.2")
AIH.pat2.2$sample <- "AIH.pat2.2"


#### 1.1.4 AIH Pat 3

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_4__20240620__110232/")
AIH.pat3 <- LoadXenium(path, fov = "fov")
AIH.pat3 <- subset(AIH.pat3, subset = nCount_Xenium > 0)
AIH.pat3 <- SCTransform(AIH.pat3, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat3),]
csv.cells <- csv.cells[match(colnames(AIH.pat3),csv.cells$cell_id),]
all(colnames(AIH.pat3)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat3 <- AddMetaData(AIH.pat3,csv.cells)
AIH.pat3.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat3",]
AIH.pat3.meta <- AIH.pat3.meta[,colnames(AIH.pat3.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat3 <- AddMetaData(AIH.pat3,AIH.pat3.meta)
AIH.pat3 <- subset(AIH.pat3,subset = celltype.simple!="other")
AIH.pat3 <- RenameCells(AIH.pat3,add.cell.id = "AIH.pat3")
AIH.pat3$sample <- "AIH.pat3"

#### 1.1.5 AIH Pat 4

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_1__20240620__110232/")
AIH.pat4 <- LoadXenium(path, fov = "fov")
AIH.pat4 <- subset(AIH.pat4, subset = nCount_Xenium > 0) 
AIH.pat4 <- SCTransform(AIH.pat4, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat4),]
csv.cells <- csv.cells[match(colnames(AIH.pat4),csv.cells$cell_id),]
all(colnames(AIH.pat4)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat4 <- AddMetaData(AIH.pat4,csv.cells)
AIH.pat4.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat4",]
AIH.pat4.meta <- AIH.pat4.meta[,colnames(AIH.pat4.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat4 <- AddMetaData(AIH.pat4,AIH.pat4.meta)
AIH.pat4 <- subset(AIH.pat4,subset = celltype.simple!="other")
AIH.pat4 <- RenameCells(AIH.pat4,add.cell.id = "AIH.pat4")
AIH.pat4$sample <- "AIH.pat4"


#### 1.1.6 AIH Pat 5

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_2__20240620__110232/")
AIH.pat5 <- LoadXenium(path, fov = "fov")
AIH.pat5 <- subset(AIH.pat5, subset = nCount_Xenium > 0)
AIH.pat5 <- SCTransform(AIH.pat5, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat5),]
csv.cells <- csv.cells[match(colnames(AIH.pat5),csv.cells$cell_id),]
all(colnames(AIH.pat5)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat5 <- AddMetaData(AIH.pat5,csv.cells)
AIH.pat5.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat5",]
AIH.pat5.meta <- AIH.pat5.meta[,colnames(AIH.pat5.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat5 <- AddMetaData(AIH.pat5,AIH.pat5.meta)
AIH.pat5 <- subset(AIH.pat5,subset = celltype.simple!="other")
AIH.pat5 <- RenameCells(AIH.pat5,add.cell.id = "AIH.pat5")
AIH.pat5$sample <- "AIH.pat5"


#### 1.1.7 AIH Pat 6
## Cut this so only the two pieces at the top are included (the bottom piece is a control biopsy)

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_3__20240620__110232/")
AIH.pat6 <- LoadXenium(path, fov = "fov")
AIH.pat6 <- subset(AIH.pat6, subset = nCount_Xenium > 0)
AIH.pat6 <- SCTransform(AIH.pat6, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(AIH.pat6),]
csv.cells <- csv.cells[match(colnames(AIH.pat6),csv.cells$cell_id),]
all(colnames(AIH.pat6)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

AIH.pat6 <- AddMetaData(AIH.pat6,csv.cells)
AIH.pat6 <- subset(AIH.pat6, subset = y_centroid<4500)
AIH.pat6.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat6",]
AIH.pat6.meta <- AIH.pat6.meta[,colnames(AIH.pat6.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
AIH.pat6 <- AddMetaData(AIH.pat6,AIH.pat6.meta)
AIH.pat6 <- subset(AIH.pat6,subset = celltype.simple!="other")
AIH.pat6 <- RenameCells(AIH.pat6,add.cell.id = "AIH.pat6")
AIH.pat6$sample <- "AIH.pat6"

### 1.2 Load all Control biopsies

#### 1.2.1 Control Pat 2

path <- paste0(working.dir,"/rawData/output-XETG00088__0021644__Region_6__20240620__110232/")
Control.pat2 <- LoadXenium(path, fov = "fov")
Control.pat2 <- subset(Control.pat2, subset = nCount_Xenium > 0)
Control.pat2 <- SCTransform(Control.pat2, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(Control.pat2),]
csv.cells <- csv.cells[match(colnames(Control.pat2),csv.cells$cell_id),]
all(colnames(Control.pat2)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

Control.pat2 <- AddMetaData(Control.pat2,csv.cells)
Control.pat2.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat2",]
Control.pat2.meta <- Control.pat2.meta[,colnames(Control.pat2.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
Control.pat2 <- AddMetaData(Control.pat2,Control.pat2.meta)
Control.pat2 <- subset(Control.pat2,subset = celltype.simple!="other")
Control.pat2 <- RenameCells(Control.pat2,add.cell.id = "Control.pat2")
Control.pat2$sample <- "Control.pat2"


#### 1.2.2 Control Pat 3
## Cut this so only the piece at the bottom is included (the top pieces are from an AIH biopsy)

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_3__20240620__110232/")
Control.pat3 <- LoadXenium(path, fov = "fov")
Control.pat3 <- subset(Control.pat3, subset = nCount_Xenium > 0)
Control.pat3 <- SCTransform(Control.pat3, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(Control.pat3),]
csv.cells <- csv.cells[match(colnames(Control.pat3),csv.cells$cell_id),]
all(colnames(Control.pat3)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

Control.pat3 <- AddMetaData(Control.pat3,csv.cells)
Control.pat3 <- subset(Control.pat3, subset = y_centroid>4500)
Control.pat3.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat3",]
Control.pat3.meta <- Control.pat3.meta[,colnames(Control.pat3.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
Control.pat3 <- AddMetaData(Control.pat3,Control.pat3.meta)
Control.pat3 <- subset(Control.pat3,subset = celltype.simple!="other")
Control.pat3 <- RenameCells(Control.pat3,add.cell.id = "Control.pat3")
Control.pat3$sample <- "Control.pat3"


#### 1.2.3 Control Pat 5

path <- paste0(working.dir,"/rawData/output-XETG00088__0021653__Region_5__20240620__110232/")
Control.pat5 <- LoadXenium(path, fov = "fov")
Control.pat5 <- subset(Control.pat5, subset = nCount_Xenium > 0)
Control.pat5 <- SCTransform(Control.pat5, assay = "Xenium")
csv.cells <- read.csv(paste0(path,"cells.csv"))
csv.cells <- csv.cells[csv.cells$cell_id %in% colnames(Control.pat5),]
csv.cells <- csv.cells[match(colnames(Control.pat5),csv.cells$cell_id),]
all(colnames(Control.pat5)==csv.cells$cell_id)
rownames(csv.cells) <- csv.cells$cell_id

Control.pat5 <- AddMetaData(Control.pat5,csv.cells)
Control.pat5.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat5",]
Control.pat5.meta <- Control.pat5.meta[,colnames(Control.pat5.meta) %in% c("celltype","celltype.TRM1","celltype.simple")]
Control.pat5 <- AddMetaData(Control.pat5,Control.pat5.meta)
Control.pat5 <- subset(Control.pat5,subset = celltype.simple!="other")
Control.pat5 <- RenameCells(Control.pat5,add.cell.id = "Control.pat5")
Control.pat5$sample <- "Control.pat5"


## 2. Merge all AIH replicates into one cell chat object

conversion.factor = 1
spot.size = 10 # use the typical human cell size

data.input1 = Seurat::GetAssayData(AIH.pat1, layer = "data", assay = "SCT") # normalized data matrix
meta1 = data.frame(labels = AIH.pat1$celltype.simple, samples = "AIH.pat1") # manually create a dataframe consisting
meta1$samples <- factor(meta1$samples)
spatial.locs1 = Seurat::GetTissueCoordinates(AIH.pat1, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs1 <- spatial.locs1[,1:2]
rownames(spatial.locs1) <- colnames(data.input1)
spatial.factors1 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input2.1 = Seurat::GetAssayData(AIH.pat2.1, layer = "data", assay = "SCT") # normalized data matrix
meta2.1 = data.frame(labels = AIH.pat2.1$celltype.simple, samples = "AIH.pat2.1") # manually create a dataframe consisting
meta2.1$samples <- factor(meta2.1$samples)
spatial.locs2.1 = Seurat::GetTissueCoordinates(AIH.pat2.1, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs2.1 <- spatial.locs2.1[,1:2]
## Add 10000 to y coordinates to prevent overlap
spatial.locs2.1$y <- spatial.locs2.1$y + 10000
rownames(spatial.locs2.1) <- colnames(data.input2.1)
spatial.factors2.1 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input2.2 = Seurat::GetAssayData(AIH.pat2.2, layer = "data", assay = "SCT") # normalized data matrix
meta2.2 = data.frame(labels = AIH.pat2.2$celltype.simple, samples = "AIH.pat2.2") # manually create a dataframe consisting
meta2.2$samples <- factor(meta2.2$samples)
spatial.locs2.2 = Seurat::GetTissueCoordinates(AIH.pat2.2, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs2.2 <- spatial.locs2.2[,1:2]
## Add 20000 to y coordinates to prevent overlap
spatial.locs2.2$y <- spatial.locs2.2$y + 20000
rownames(spatial.locs2.2) <- colnames(data.input2.2)
spatial.factors2.2 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input3 = Seurat::GetAssayData(AIH.pat3, layer = "data", assay = "SCT") # normalized data matrix
meta3 = data.frame(labels = AIH.pat3$celltype.simple, samples = "AIH.pat3") # manually create a dataframe consisting
meta3$samples <- factor(meta3$samples)
spatial.locs3 = Seurat::GetTissueCoordinates(AIH.pat3, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs3 <- spatial.locs3[,1:2]
## Add 30000 to y coordinates to prevent overlap
spatial.locs3$y <- spatial.locs3$y + 30000
rownames(spatial.locs3) <- colnames(data.input3)
spatial.factors3 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input4 = Seurat::GetAssayData(AIH.pat4, layer = "data", assay = "SCT") # normalized data matrix
meta4 = data.frame(labels = AIH.pat4$celltype.simple, samples = "AIH.pat4") # manually create a dataframe consisting
meta4$samples <- factor(meta4$samples)
spatial.locs4 = Seurat::GetTissueCoordinates(AIH.pat4, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs4 <- spatial.locs4[,1:2]
## Add 40000 to y coordinates to prevent overlap
spatial.locs4$y <- spatial.locs4$y + 40000
rownames(spatial.locs4) <- colnames(data.input4)
spatial.factors4 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input5 = Seurat::GetAssayData(AIH.pat5, layer = "data", assay = "SCT") # normalized data matrix
meta5 = data.frame(labels = AIH.pat5$celltype.simple, samples = "AIH.pat5") # manually create a dataframe consisting
meta5$samples <- factor(meta5$samples)
spatial.locs5 = Seurat::GetTissueCoordinates(AIH.pat5, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs5 <- spatial.locs5[,1:2]
## Add 50000 to y coordinates to prevent overlap
spatial.locs5$y <- spatial.locs5$y + 50000
rownames(spatial.locs5) <- colnames(data.input5)
spatial.factors5 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input6 = Seurat::GetAssayData(AIH.pat6, layer = "data", assay = "SCT") # normalized data matrix
meta6 = data.frame(labels = AIH.pat6$celltype.simple, samples = "AIH.pat6") # manually create a dataframe consisting
meta6$samples <- factor(meta6$samples)
spatial.locs6 = Seurat::GetTissueCoordinates(AIH.pat6, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs6 <- spatial.locs6[,1:2]
## Add 60000 to y coordinates to prevent overlap
spatial.locs6$y <- spatial.locs6$y + 60000
rownames(spatial.locs6) <- colnames(data.input6)
spatial.factors6 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input.AIH <- do.call(cbind,list(data.input1,data.input2.1,data.input2.2,data.input3,data.input4,data.input5,data.input6))
meta.AIH <- do.call(rbind,list(meta1,meta2.1,meta2.2,meta3,meta4,meta5,meta6))
meta.AIH$labels <- droplevels(meta.AIH$labels, exclude = setdiff(levels(meta.AIH$labels),unique(meta.AIH$labels)))
spatial.locs.AIH <- do.call(rbind,list(spatial.locs1,spatial.locs2.1,spatial.locs2.2,spatial.locs3,spatial.locs4,spatial.locs5,spatial.locs6))
spatial.factors.AIH <- do.call(rbind,list(spatial.factors1,spatial.factors2.1,spatial.factors2.2,spatial.factors3,spatial.factors4,spatial.factors5,spatial.factors6))
cellchat.AIH.subset <- createCellChat(object = data.input.AIH, meta = meta.AIH, group.by = "labels",datatype = "spatial", coordinates = spatial.locs.AIH, spatial.factors = spatial.factors.AIH)
cellchat.AIH.subset


## 3. Merge all Control replicates

## Merge all Control replicates into one cell chat object

conversion.factor = 1
spot.size = 10 # use the typical human cell size                                                     

data.input2 = Seurat::GetAssayData(Control.pat2, layer = "data", assay = "SCT") # normalized data matrix
meta2 = data.frame(labels = Control.pat2$celltype.simple, samples = "Control.pat2") # manually create a dataframe consisting
meta2$samples <- factor(meta2$samples)
spatial.locs2 = Seurat::GetTissueCoordinates(Control.pat2, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs2 <- spatial.locs2[,1:2]
## Add 20000 to y coordinates to prevent overlap
spatial.locs2$y <- spatial.locs2$y + 10000
spatial.locs2$x <- spatial.locs2$x + 10000
rownames(spatial.locs2) <- colnames(data.input2)
spatial.factors2 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input3 = Seurat::GetAssayData(Control.pat3, layer = "data", assay = "SCT") # normalized data matrix
meta3 = data.frame(labels = Control.pat3$celltype.simple, samples = "Control.pat3") # manually create a dataframe consisting
meta3$samples <- factor(meta3$samples)
spatial.locs3 = Seurat::GetTissueCoordinates(Control.pat3, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs3 <- spatial.locs3[,1:2]
## Add 30000 to y coordinates to prevent overlap
spatial.locs3$y <- spatial.locs3$y + 20000
spatial.locs3$x <- spatial.locs3$x + 20000
rownames(spatial.locs3) <- colnames(data.input3)
spatial.factors3 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input5 = Seurat::GetAssayData(Control.pat5, layer = "data", assay = "SCT") # normalized data matrix
meta5 = data.frame(labels = Control.pat5$celltype.simple, samples = "Control.pat5") # manually create a dataframe consisting
meta5$samples <- factor(meta5$samples)
spatial.locs5 = Seurat::GetTissueCoordinates(Control.pat5, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs5 <- spatial.locs5[,1:2]
## Add 50000 to y coordinates to prevent overlap
spatial.locs5$y <- spatial.locs5$y + 40000
spatial.locs5$x <- spatial.locs5$x + 10000
rownames(spatial.locs5) <- colnames(data.input5)
spatial.factors5 = data.frame(ratio = conversion.factor, tol = spot.size/2)

data.input.Control <- do.call(cbind,list(data.input2,data.input3,data.input5))
meta.Control <- do.call(rbind,list(meta2,meta3,meta5))
meta.Control$labels <- droplevels(meta.Control$labels, exclude = setdiff(levels(meta.Control$labels),unique(meta.Control$labels)))
spatial.locs.Control <- do.call(rbind,list(spatial.locs2,spatial.locs3,spatial.locs5))
spatial.factors.Control <- do.call(rbind,list(spatial.factors2,spatial.factors3,spatial.factors5))

cellchat.Control.subset <- createCellChat(object = data.input.Control, meta = meta.Control, group.by = "labels",datatype = "spatial", coordinates = spatial.locs.Control, spatial.factors = spatial.factors.Control)
cellchat.Control.subset

saveRDS(cellchat.AIH.subset,paste0(working.dir,"/SeuratObjects/AIH_CellChat.rds"))
saveRDS(cellchat.Control.subset,paste0(working.dir,"/SeuratObjects/Control_CellChat.rds"))

