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
Control.combined.integrated <- readRDS(file = paste0(working.dir,"/SeuratObjects/Control_integrated.rds"))


## Extended Figure 6B - AIH UMAP -----------------------------------------------

pdf(paste0(working.dir,"/results/AIH_integrated_UMAP_renamed.pdf"),width = 16,height = 10)
DimPlot(AIH.combined.integrated,label = TRUE) & coord_fixed() & NoAxes()
dev.off()

## Extended Figure 6B - AIH Dotplot --------------------------------------------

genes.heatmap <- c("CYP2E1","GLUL","SDS","HAL","CXCL10","IL32","KRT7","EPCAM","FHL2","PDGFRA","DNASE1L3","ADGRL4","LYVE1","MMRN1","MYH11","MYLK","CSF1R","MPEG1","CD5L","MARCO","CPA3","CTSG","IDO1","XCR1","LILRA4","KLRC1","KLRB1","CD3E","CD8A","CD4","CD79A","CD19","SLAMF7","MZB1","MKI67","TOP2A")

pdf(paste0(working.dir,"/results/AIH_integrated_cluster_markers_Dotplot.pdf"),width = 16,height = 10)
DotPlot(AIH.combined.integrated,genes.heatmap,group.by = "celltype") + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) + coord_fixed()
dev.off()

## Extended Figure 6C - Control UMAP -------------------------------------------

pdf(paste0(working.dir,"/results/Control_integrated_UMAP_renamed.pdf"),width = 16,height = 10)
DimPlot(Control.combined.integrated,label = TRUE) & coord_fixed() & NoAxes()
dev.off()

## Extended Figure 6C - Control Dotplot ----------------------------------------

genes.heatmap <- c("CYP2E1","GLUL","SDS","HAL","CXCL10","IL32","KRT7","EPCAM","FHL2","PDGFRA","DNASE1L3","ADGRL4","LYVE1","MMRN1","MYH11","MYLK","CSF1R","MPEG1","CD5L","MARCO","KLRC1","KLRB1","CD3E","CD8A","CD4","CD79A","CD19","SLAMF7","MZB1","MKI67","TOP2A")

pdf(paste0(working.dir,"/results/Control_integrated_cluster_markers_Dotplot.pdf"),width = 16,height = 10)
DotPlot(Control.combined.integrated,genes.heatmap,group.by = "celltype") + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) + coord_fixed()
dev.off()

## Extended Figure 6D - Proportions --------------------------------------------

AIH.combined.integrated$celltype.simplified <- as.character(AIH.combined.integrated$celltype.TRM1)
AIH.combined.integrated$celltype.simplified[grepl("Hepat",AIH.combined.integrated$celltype.TRM1)] <- "Hepatocytes"
AIH.combined.integrated$celltype.simplified[grepl("TEM",AIH.combined.integrated$celltype.TRM1)] <- "Tcells.CD4.TEM"
AIH.combined.integrated$celltype.simplified[grepl("aaCD8",AIH.combined.integrated$celltype.TRM1)] <- "Tcells.aaCD8"
AIH.combined.integrated$celltype.simplified[grepl("TRM1",AIH.combined.integrated$celltype.TRM1)] <- "Tcells.CD4.TRM1"
AIH.combined.integrated$celltype.simplified[grepl("EC.",AIH.combined.integrated$celltype.TRM1)] <- "EC"
AIH.combined.integrated$celltype.simplified[AIH.combined.integrated$celltype %in% c("Kupffer.cells","Mac.Monocytes","cDC1s","pDCs","Granulocytes")] <- "Myeloid.cells"
AIH.combined.integrated$celltype.simplified[AIH.combined.integrated$celltype.TRM1 %in% c("dividing.cells","undefined","SMC","Neutrophils","ILC","NK.cells","Tcells.CD4")] <- "other"

Control.combined.integrated$celltype.simplified <- as.character(Control.combined.integrated$celltype.TRM1)
Control.combined.integrated$celltype.simplified[grepl("Hepat",Control.combined.integrated$celltype.TRM1)] <- "Hepatocytes"
Control.combined.integrated$celltype.simplified[grepl("TEM",Control.combined.integrated$celltype.TRM1)] <- "Tcells.CD4.TEM"
Control.combined.integrated$celltype.simplified[grepl("aaCD8",Control.combined.integrated$celltype.TRM1)] <- "Tcells.aaCD8"
Control.combined.integrated$celltype.simplified[grepl("TRM1",Control.combined.integrated$celltype.TRM1)] <- "Tcells.CD4.TRM1"
Control.combined.integrated$celltype.simplified[grepl("EC.",Control.combined.integrated$celltype.TRM1)] <- "EC"
Control.combined.integrated$celltype.simplified[Control.combined.integrated$celltype %in% c("Kupffer.cells","Mac.Monocytes","cDC1s","pDCs","Granulocytes")] <- "Myeloid.cells"
Control.combined.integrated$celltype.simplified[Control.combined.integrated$celltype.TRM1 %in% c("dividing.cells","undefined","SMC","Neutrophils","ILC","NK.cells","Tcells.CD4")] <- "other"

prop.celltypes.AIH <- as.data.frame(prop.table(table(AIH.combined.integrated$celltype.simplified,AIH.combined.integrated$patient),margin = 2))
colnames(prop.celltypes.AIH) <- c("celltype.simplified","patient","proportion")

prop.celltypes.Control <- as.data.frame(prop.table(table(Control.combined.integrated$celltype.simplified,Control.combined.integrated$patient),margin = 2))
colnames(prop.celltypes.Control) <- c("celltype.simplified","patient","proportion")

prop.celltypes.all <- rbind(prop.celltypes.AIH,prop.celltypes.Control)

prop.celltypes.all$disease <- ifelse(grepl("AIH",prop.celltypes.all$patient),"AIH","Control")
colnames(prop.celltypes.all) <- c("celltype.simplified","patient","proportion","disease")

prop.celltypes.all$celltype.simplified <- factor(prop.celltypes.all$celltype.simplified,levels = c("Tcells.CD4.TRM1","Tcells.aaCD8","Tcells.CD4.TEM","Tcells.CD8","Myeloid.cells","B.cells","Plasma.cells","Hepatocytes","Cholangiocytes","EC","Fibroblasts","other"))
prop.celltypes.all$percentage <- prop.celltypes.all$proportion*100

pdf(paste0(working.dir,"/results/Proportions_celltypes.pdf"),width = 20,height = 3.5)
stat.test <- prop.celltypes.all %>% group_by(celltype.simplified) %>% t_test(percentage ~ disease) %>% adjust_pvalue(method = "bonferroni")
p1 <- ggbarplot(prop.celltypes.all, x = "disease", y = "percentage", color = "disease",add = c("mean_se","jitter"), facet.by = "celltype.simplified",scales = "free",nrow = 1) +
  scale_color_manual(values=c("darkred","darkblue")) +
  theme_classic()  +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.07)))
stat.test <- stat.test %>% add_xy_position(fun = "max", x = "disease")
p1 <- p1 + stat_pvalue_manual(stat.test)
print(p1)
dev.off()


pdf(paste0(working.dir,"/results/Proportions_celltypes_no_pval.pdf"),width = 20,height = 3.5)
p1 <- ggbarplot(prop.celltypes.all, x = "disease", y = "percentage", color = "disease",add = c("mean_se","jitter"), facet.by = "celltype.simplified",scales = "free",nrow = 1) +
  scale_color_manual(values=c("darkred","darkblue")) +
  theme_classic()
print(p1)
dev.off()

## Extended Figure 6E - Hepatocytes --------------------------------------------

## add a patient column to make sure that the two biopsies from pat2 are counted together
AIH.combined.integrated$patient <- as.character(AIH.combined.integrated$sample)
AIH.combined.integrated$patient[grepl("AIH.pat2",AIH.combined.integrated$sample)] <- "AIH.pat2"

AIH.hepat.all <- as.data.frame(prop.table(table(AIH.combined.integrated$celltype[grepl("Hepat",AIH.combined.integrated$celltype)],AIH.combined.integrated$patient[grepl("Hepat",AIH.combined.integrated$celltype)]),margin = 2))
AIH.hepat.all <- AIH.hepat.all[AIH.hepat.all$Freq>0,]
colnames(AIH.hepat.all) <- c("celltype","patient","proportion")

Control.combined.integrated$patient <- Control.combined.integrated$sample
Control.hepat.all <- as.data.frame(prop.table(table(Control.combined.integrated$celltype[grepl("Hepat",Control.combined.integrated$celltype)],Control.combined.integrated$patient[grepl("Hepat",Control.combined.integrated$celltype)]),margin = 2))
Control.hepat.all <- Control.hepat.all[Control.hepat.all$Freq>0,]
colnames(Control.hepat.all) <- c("celltype","patient","proportion")

Hepat.all <- rbind(AIH.hepat.all,Control.hepat.all)
Hepat.all$patient <- factor(Hepat.all$patient)
Hepat.all$celltype <- factor(Hepat.all$celltype,levels=c("Hepatocytes.central","Hepatocytes.midzone.1","Hepatocytes.midzone.2","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory"))

pdf(paste0(working.dir,"/results/Hepatocyte_distribution.pdf"),width = 6,height = 4)
ggplot(Hepat.all,aes(x=patient,y=proportion,fill=celltype)) +
  geom_col() +
  scale_fill_manual(values = rev(c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5"))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

## Extended Figure 6F - Neighborhood analysis ----------------------------------

## Rationale: (1) split objects into individual biopsies. Then, for each biopsy: 
## (1.1) Iterate over all cells and find the 5 closest cells for each using knn. Note: use k=6, since the result will include the cell itself. Then exclude that cell from the results and only ## keep the 5 closest cells. Store the celltypes of these 5 cells in a list at the same index as that cell in the biopsy meta data.
## (1.2) Iterate over all celltypes to create a summary. For each celltype, find the indeces of these cells in the biopsy meta data. Then extract the neighboring cells of only those indeces from the previous list and tally them.

## Control biopsies
Control.combined.integrated$celltype.simplified <- factor(Control.combined.integrated$celltype.simplified,levels = c("Hepatocytes","Cholangiocytes","Fibroblasts","EC","Myeloid.cells","Tcells.CD4.TEM","Tcells.CD4.TRM1","Tcells.CD8","Tcells.aaCD8","B.cells","Plasma.cells","other"))

Control.pat2.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat2",]
Control.pat2.neighborhoods <- list()
for (i in 1:nrow(Control.pat2.meta)){
  knn_result <- get.knnx(data.frame(x=Control.pat2.meta$x_centroid, y=Control.pat2.meta$y_centroid), data.frame(x=Control.pat2.meta$x_centroid, y=Control.pat2.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  Control.pat2.neighborhoods[[i]] <- list(celltypes = Control.pat2.meta$celltype.simplified[knn.numbers])
}
Control.pat2.summary <-list()
for (i in 1:nlevels(Control.pat2.meta$celltype.simplified)){
  indices.celltype.simplified <- which(Control.pat2.meta$celltype.simplified==levels(Control.pat2.meta$celltype.simplified)[i])
  Control.pat2.summary[[i]] <-
    list(unlist(Control.pat2.neighborhoods[indices.celltype.simplified]))
}

Control.pat3.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat3",]
Control.pat3.neighborhoods <- list()
for (i in 1:nrow(Control.pat3.meta)){
  knn_result <- get.knnx(data.frame(x=Control.pat3.meta$x_centroid, y=Control.pat3.meta$y_centroid), data.frame(x=Control.pat3.meta$x_centroid, y=Control.pat3.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  Control.pat3.neighborhoods[[i]] <- list(celltypes = Control.pat3.meta$celltype.simplified[knn.numbers])
}
Control.pat3.summary <-list()
for (i in 1:nlevels(Control.pat3.meta$celltype.simplified)){
  indices.celltype.simplified <- which(Control.pat3.meta$celltype.simplified==levels(Control.pat3.meta$celltype.simplified)[i])
  Control.pat3.summary[[i]] <-
    list(unlist(Control.pat3.neighborhoods[indices.celltype.simplified]))
}

Control.pat5.meta <- Control.combined.integrated@meta.data[Control.combined.integrated$sample=="Control.pat5",]
Control.pat5.neighborhoods <- list()
for (i in 1:nrow(Control.pat5.meta)){
  knn_result <- get.knnx(data.frame(x=Control.pat5.meta$x_centroid, y=Control.pat5.meta$y_centroid), data.frame(x=Control.pat5.meta$x_centroid, y=Control.pat5.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  Control.pat5.neighborhoods[[i]] <- list(celltypes = Control.pat5.meta$celltype.simplified[knn.numbers])
}
Control.pat5.summary <-list()
for (i in 1:nlevels(Control.pat5.meta$celltype.simplified)){
  indices.celltype.simplified <- which(Control.pat5.meta$celltype.simplified==levels(Control.pat5.meta$celltype.simplified)[i])
  Control.pat5.summary[[i]] <-
    list(unlist(Control.pat5.neighborhoods[indices.celltype.simplified]))
}

# Summarize Control biopsies
Control.prop.knn <- data.frame(celltype.simplified = character(0),neighbor.celltype.simplified = character(0),sample = character(0),proportion = numeric(0))
for (i in 1:nlevels(Control.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(Control.pat2.summary[[i]])))
  Control.prop.knn <- rbind(Control.prop.knn,data.frame(celltype.simplified = rep(levels(Control.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("Control.pat2",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
for (i in 1:nlevels(Control.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(Control.pat3.summary[[i]])))
  Control.prop.knn <- rbind(Control.prop.knn,data.frame(celltype.simplified = rep(levels(Control.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("Control.pat3",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
for (i in 1:nlevels(Control.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(Control.pat5.summary[[i]])))
  Control.prop.knn <- rbind(Control.prop.knn,data.frame(celltype.simplified = rep(levels(Control.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("Control.pat5",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}

## AIH biopsies

AIH.combined.integrated$celltype.simplified <- factor(AIH.combined.integrated$celltype.simplified,levels = c("Hepatocytes","Cholangiocytes","Fibroblasts","EC","Myeloid.cells","Tcells.CD4.TEM","Tcells.CD4.TRM1","Tcells.CD8","Tcells.aaCD8","B.cells","Plasma.cells","other"))

AIH.pat1.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat1",]
AIH.pat1.neighborhoods <- list()
for (i in 1:nrow(AIH.pat1.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat1.meta$x_centroid, y=AIH.pat1.meta$y_centroid), data.frame(x=AIH.pat1.meta$x_centroid, y=AIH.pat1.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat1.neighborhoods[[i]] <- list(celltypes = AIH.pat1.meta$celltype.simplified[knn.numbers])
}
AIH.pat1.summary <-list()
for (i in 1:nlevels(AIH.pat1.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat1.meta$celltype.simplified==levels(AIH.pat1.meta$celltype.simplified)[i])
  AIH.pat1.summary[[i]] <-
    list(unlist(AIH.pat1.neighborhoods[indices.celltype.simplified]))
}

AIH.pat2.1.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.1",]
AIH.pat2.1.neighborhoods <- list()
for (i in 1:nrow(AIH.pat2.1.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat2.1.meta$x_centroid, y=AIH.pat2.1.meta$y_centroid), data.frame(x=AIH.pat2.1.meta$x_centroid, y=AIH.pat2.1.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat2.1.neighborhoods[[i]] <- list(celltypes = AIH.pat2.1.meta$celltype.simplified[knn.numbers])
}
AIH.pat2.1.summary <-list()
for (i in 1:nlevels(AIH.pat2.1.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat2.1.meta$celltype.simplified==levels(AIH.pat2.1.meta$celltype.simplified)[i])
  AIH.pat2.1.summary[[i]] <-
    list(unlist(AIH.pat2.1.neighborhoods[indices.celltype.simplified]))
}

AIH.pat2.2.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.2",]
AIH.pat2.2.neighborhoods <- list()
for (i in 1:nrow(AIH.pat2.2.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat2.2.meta$x_centroid, y=AIH.pat2.2.meta$y_centroid), data.frame(x=AIH.pat2.2.meta$x_centroid, y=AIH.pat2.2.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat2.2.neighborhoods[[i]] <- list(celltypes = AIH.pat2.2.meta$celltype.simplified[knn.numbers])
}
AIH.pat2.2.summary <-list()
for (i in 1:nlevels(AIH.pat2.2.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat2.2.meta$celltype.simplified==levels(AIH.pat2.2.meta$celltype.simplified)[i])
  AIH.pat2.2.summary[[i]] <-
    list(unlist(AIH.pat2.2.neighborhoods[indices.celltype.simplified]))
}

AIH.pat3.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat3",]
AIH.pat3.neighborhoods <- list()
for (i in 1:nrow(AIH.pat3.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat3.meta$x_centroid, y=AIH.pat3.meta$y_centroid), data.frame(x=AIH.pat3.meta$x_centroid, y=AIH.pat3.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat3.neighborhoods[[i]] <- list(celltypes = AIH.pat3.meta$celltype.simplified[knn.numbers])
}
AIH.pat3.summary <-list()
for (i in 1:nlevels(AIH.pat3.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat3.meta$celltype.simplified==levels(AIH.pat3.meta$celltype.simplified)[i])
  AIH.pat3.summary[[i]] <-
    list(unlist(AIH.pat3.neighborhoods[indices.celltype.simplified]))
}

AIH.pat4.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat4",]
AIH.pat4.neighborhoods <- list()
for (i in 1:nrow(AIH.pat4.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat4.meta$x_centroid, y=AIH.pat4.meta$y_centroid), data.frame(x=AIH.pat4.meta$x_centroid, y=AIH.pat4.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat4.neighborhoods[[i]] <- list(celltypes = AIH.pat4.meta$celltype.simplified[knn.numbers])
}
AIH.pat4.summary <-list()
for (i in 1:nlevels(AIH.pat4.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat4.meta$celltype.simplified==levels(AIH.pat4.meta$celltype.simplified)[i])
  AIH.pat4.summary[[i]] <-
    list(unlist(AIH.pat4.neighborhoods[indices.celltype.simplified]))
}

AIH.pat5.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat5",]
AIH.pat5.neighborhoods <- list()
for (i in 1:nrow(AIH.pat5.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat5.meta$x_centroid, y=AIH.pat5.meta$y_centroid), data.frame(x=AIH.pat5.meta$x_centroid, y=AIH.pat5.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat5.neighborhoods[[i]] <- list(celltypes = AIH.pat5.meta$celltype.simplified[knn.numbers])
}
AIH.pat5.summary <-list()
for (i in 1:nlevels(AIH.pat5.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat5.meta$celltype.simplified==levels(AIH.pat5.meta$celltype.simplified)[i])
  AIH.pat5.summary[[i]] <-
    list(unlist(AIH.pat5.neighborhoods[indices.celltype.simplified]))
}

AIH.pat6.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat6",]
AIH.pat6.neighborhoods <- list()
for (i in 1:nrow(AIH.pat6.meta)){
  knn_result <- get.knnx(data.frame(x=AIH.pat6.meta$x_centroid, y=AIH.pat6.meta$y_centroid), data.frame(x=AIH.pat6.meta$x_centroid, y=AIH.pat6.meta$y_centroid)[i,], k=6)
  knn.numbers <- knn_result$nn.index[knn_result$nn.index !=i]
  AIH.pat6.neighborhoods[[i]] <- list(celltypes = AIH.pat6.meta$celltype.simplified[knn.numbers])
}
AIH.pat6.summary <-list()
for (i in 1:nlevels(AIH.pat6.meta$celltype.simplified)){
  indices.celltype.simplified <- which(AIH.pat6.meta$celltype.simplified==levels(AIH.pat6.meta$celltype.simplified)[i])
  AIH.pat6.summary[[i]] <-
    list(unlist(AIH.pat6.neighborhoods[indices.celltype.simplified]))
}

## Summarize AIH biopsies

AIH.prop.knn <- data.frame(celltype.simplified = character(0),neighbor.celltype.simplified = character(0),sample = character(0),proportion = numeric(0))
for (i in 1:nlevels(AIH.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(AIH.pat1.summary[[i]])))
  AIH.prop.knn <- rbind(AIH.prop.knn,data.frame(celltype.simplified = rep(levels(AIH.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("AIH.pat1",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
## concatenate the results for the two AIH pat2 biopsies
for (i in 1:nlevels(AIH.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(unlist(c(AIH.pat2.1.summary[[i]],AIH.pat2.2.summary[[i]])))))
  AIH.prop.knn <- rbind(AIH.prop.knn,data.frame(celltype.simplified = rep(levels(AIH.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("AIH.pat2",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
for (i in 1:nlevels(AIH.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(AIH.pat3.summary[[i]])))
  AIH.prop.knn <- rbind(AIH.prop.knn,data.frame(celltype.simplified = rep(levels(AIH.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("AIH.pat3",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
for (i in 1:nlevels(AIH.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(AIH.pat4.summary[[i]])))
  AIH.prop.knn <- rbind(AIH.prop.knn,data.frame(celltype.simplified = rep(levels(AIH.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("AIH.pat4",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
for (i in 1:nlevels(AIH.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(AIH.pat5.summary[[i]])))
  AIH.prop.knn <- rbind(AIH.prop.knn,data.frame(celltype.simplified = rep(levels(AIH.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("AIH.pat5",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
for (i in 1:nlevels(AIH.combined.integrated$celltype.simplified)){
  celltypes.i <- as.data.frame(prop.table(table(AIH.pat6.summary[[i]])))
  AIH.prop.knn <- rbind(AIH.prop.knn,data.frame(celltype.simplified = rep(levels(AIH.combined.integrated$celltype.simplified)[i],nrow(celltypes.i)),neighbor.celltype.simplified = celltypes.i$Var1,sample = rep("AIH.pat6",nrow(celltypes.i)),proportion = celltypes.i$Freq))
}
## Summarize and Plot

knn.prop.all <- rbind(AIH.prop.knn,Control.prop.knn)
knn.prop.all$disease <- ifelse(grepl("AIH",knn.prop.all$sample),"AIH","Control")
knn.prop.all$percentage <- knn.prop.all$proportion*100

pdf(paste0(working.dir,"/results/Neighboring_cells_individual.pdf"),width = 18,height = 3.5)
for (i in 1:length(unique(knn.prop.all$celltype.simplified))){
  source.celltype <- unique(knn.prop.all$celltype.simplified)[i]
  prop.red <- knn.prop.all[knn.prop.all$celltype.simplified==source.celltype,]
  stat.test <- prop.red %>% group_by(neighbor.celltype.simplified) %>% t_test(percentage ~ disease) %>% adjust_pvalue(method = "bonferroni")
  p1 <- ggbarplot(prop.red, x = "disease", y = "percentage", color = "disease",add = c("mean_se","jitter"), facet.by = "neighbor.celltype.simplified",scales = "free",nrow = 1) +
    scale_color_manual(values=c("darkred","darkblue")) +
    theme_classic()  +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(title = paste0("Neighboring cells of ",source.celltype))
  stat.test <- stat.test %>% add_xy_position(fun = "max", x = "disease")
  p1 <- p1 + stat_pvalue_manual(stat.test)
  print(p1)
}
dev.off()


pdf(paste0(working.dir,"/results/Neighboring_cells_individual_no_pval.pdf"),width = 18,height = 3.5)
for (i in 1:length(unique(knn.prop.all$celltype.simplified))){
  source.celltype <- unique(knn.prop.all$celltype.simplified)[i]
  prop.red <- knn.prop.all[knn.prop.all$celltype.simplified==source.celltype,]
  # stat.test <- prop.red %>% group_by(neighbor.celltype.simplified) %>% t_test(percentage ~ disease) %>% adjust_pvalue(method = "bonferroni")
  p1 <- ggbarplot(prop.red, x = "disease", y = "percentage", color = "disease",add = c("mean_se","jitter"), facet.by = "neighbor.celltype.simplified",scales = "free",nrow = 1) +
    scale_color_manual(values=c("darkred","darkblue")) +
    theme_classic()  +
    # scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(title = paste0("Neighboring cells of ",source.celltype))
  # stat.test <- stat.test %>% add_xy_position(fun = "max", x = "disease")
  # p1 <- p1 + stat_pvalue_manual(stat.test)
  print(p1)
}
dev.off()


