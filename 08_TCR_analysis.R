## TCR analysis ##

library(Seurat)
library(dplyr)
library(ggplot2)
library(scRepertoire)


PBMC_09 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_PBMC_AIH09.rds")
PBMC_09$protocol <- "SCS"

PBMC_13 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_PBMC_AIH13.rds")
PBMC_13$protocol <- "SCS"

PBMC_15 <- readRDS("./Single_sample_preprocess/RDS/prefiltered_PBMC_AIH15.rds")
PBMC_15$protocol <- "SCS"

PBMC <- merge(PBMC_09, y = c(PBMC_13, 
                              PBMC_15), 
               add.cell.ids = c("PBMC_AIH09", 
                                "PBMC_AIH13", 
                                "PBMC_AIH15"))
PBMC.AIH09.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/PBMC_AIH09_singlet_barcodes.txt",header = FALSE)$V1
PBMC.AIH13.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/PBMC_AIH13_singlet_barcodes.txt",header = FALSE)$V1
PBMC.AIH15.singlet.barcodes <- read.delim("./Single_sample_preprocess/Doublet_finder_table/PBMC_AIH15_singlet_barcodes.txt",header = FALSE)$V1

all.SCS.singlet.barcodes <- c(PBMC.AIH09.singlet.barcodes, 
                              PBMC.AIH13.singlet.barcodes,
                              PBMC.AIH15.singlet.barcodes)

PBMC$DF.status <- ifelse(rownames(PBMC@meta.data) %in% all.SCS.singlet.barcodes,"Singlet","Doublet")

PBMC_singlets <- subset(PBMC, subset = DF.status == "Singlet")
PBMC_singlets$tissue <- "PBMC"
PBMC_singlets$DF.status <- NULL


#Import liver T cells from AIH atlas
AIH.combined <- readRDS(file='./Integration/RDS/3_SCS_SNS_atlas.rds')
Liver_Tcell <- subset(x = AIH.combined, idents = c("CD4", "CD8"))

Liver_Tcell$annotation <- NULL
Liver_Tcell$tissue <- "Liver"

Liver_Tcell_SCS <- subset(Liver_Tcell, subset = protocol == "SCS")

Tcells.all <- merge(PBMC_singlets, y = Liver_Tcell_SCS)

sample.list <- SplitObject(Tcells.all, split.by = "orig.ident")

sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)

Tcells.all.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)

Tcells.all.combined <- IntegrateData(anchorset = Tcells.all.anchors)

DefaultAssay(Tcells.all.combined) <- "integrated"
Tcells.all.combined <- ScaleData(Tcells.all.combined, verbose = FALSE)
Tcells.all.combined <- RunPCA(Tcells.all.combined, npcs = 30, verbose = FALSE)
Tcells.all.combined <- RunUMAP(Tcells.all.combined, reduction = "pca", dims = 1:30)
Tcells.all.combined <- FindNeighbors(Tcells.all.combined, reduction = "pca", dims = 1:30)
Tcells.all.combined <- FindClusters(Tcells.all.combined, resolution = 0.6)

DefaultAssay(Tcells.all.combined) <- "RNA"
Tcells.all.combined <- NormalizeData(Tcells.all.combined, normalization.method = "LogNormalize", margin = 1, assay = "RNA")
Tcells.all.combined <- ScaleData(Tcells.all.combined) 

DefaultAssay(Tcells.all.combined) <- "CITE"
Tcells.all.combined <- NormalizeData(Tcells.all.combined, normalization.method = "CLR", margin = 2, assay = "CITE")
Tcells.all.combined <- ScaleData(Tcells.all.combined) 

## Subset to patients with matching blood and liver samples
Tcells.all.combined.subset <- subset(Tcells.all.combined, idents= c("PBMC-A09", "PBMC-A13", "PBMC-A15", "A09", "A13", "A15"))

s1 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH9_TCR/filtered_contig_annotations.csv")
s2 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH13_TCR/filtered_contig_annotations.csv")
s3 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH15_TCR/filtered_contig_annotations.csv")
s4 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH9_PBMC_TCR/filtered_contig_annotations.csv")
s5 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH13_PBMC_TCR/filtered_contig_annotations.csv")
s6 <- read.csv ("./Single_sample_Cellranger/SC-Seq-Cellranger/AIH15_PBMC_TCR/filtered_contig_annotations.csv")

ob.list <- list ()
samplelist <- c("AIH09", "AIH13", "AIH15", "PBMC_AIH09", "PBMC_AIH13", "PBMC_AIH15")
ob.list[[1]] <- s1
ob.list[[2]] <- s2
ob.list[[3]] <- s3
ob.list[[4]] <- s4
ob.list[[5]] <- s5
ob.list[[6]] <- s6

combined <- combineTCR(ob.list,
                       sample=samplelist,
                       filterMulti=FALSE)

Tcells.all.combined.subset <- combineExpression(combined, 
                         Tcells.all.combined.subset,
                         cloneCall = "gene+nt",
                         proportion = F,
                         cloneType = c(None = 0, Single =1, Small = 10, Medium = 70, Large =200),
                         filterNA= T)

saveRDS(Tcells.all.combined.subset,file = "./Integration/RDS/9_PBMC_liver_Tcell_clustered.rds")