## Figure 4 ##

library(Seurat) # packageVersion("Seurat") # 5.0.3
library(tidyverse) # packageVersion("tidyverse") # 2.0.0
library(org.Hs.eg.db) # packageVersion("org.Hs.eg.db") # 3.16.0
library(KEGGREST) # packageVersion("KEGGREST") # 1.38.0
library(ggrepel) # packageVersion("ggrepel") # 0.9.5

## Insert directory here which contains the following folders:
## rawData #(contains Xenium raw data)
## results
## doc
## SeuratObjects
# working.dir <- "path/to/working/directory"

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


AIH.combined.integrated <- readRDS(file = paste0(working.dir,"/SeuratObjects/AIH_integrated.rds"))


## Figure 5H -------------------------------------------------------------------

## 1. Load AIH objects and determine distance to closest TNF producer for each sample

AIH.combined.integrated$TNF.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="TNF",]
AIH.combined.integrated$TNF.status <- ifelse(AIH.combined.integrated$TNF.expr>0,"TNF.pos","TNF.neg")

## AIH.pat1

AIH.pat1.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat1",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat1.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat1",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat1.meta)){
  x.pos.cell <- AIH.pat1.meta$x_centroid[i]
  y.pos.cell <- AIH.pat1.meta$y_centroid[i]
  annot.df <- AIH.pat1.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat1.distances <- data.frame(cell_id=AIH.pat1.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat1.distances) <- paste0("AIH.pat1_",AIH.pat1.distances$cell_id)

## AIH.pat2.1

AIH.pat2.1.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat2.1",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat2.1.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.1",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat2.1.meta)){
  x.pos.cell <- AIH.pat2.1.meta$x_centroid[i]
  y.pos.cell <- AIH.pat2.1.meta$y_centroid[i]
  annot.df <- AIH.pat2.1.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat2.1.distances <- data.frame(cell_id=AIH.pat2.1.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat2.1.distances) <- paste0("AIH.pat2.1_",AIH.pat2.1.distances$cell_id)


## AIH.pat2.1

AIH.pat2.2.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat2.2",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat2.2.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat2.2",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat2.2.meta)){
  x.pos.cell <- AIH.pat2.2.meta$x_centroid[i]
  y.pos.cell <- AIH.pat2.2.meta$y_centroid[i]
  annot.df <- AIH.pat2.2.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat2.2.distances <- data.frame(cell_id=AIH.pat2.2.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat2.2.distances) <- paste0("AIH.pat2.2_",AIH.pat2.2.distances$cell_id)


## AIH.pat3

AIH.pat3.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat3",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat3.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat3",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat3.meta)){
  x.pos.cell <- AIH.pat3.meta$x_centroid[i]
  y.pos.cell <- AIH.pat3.meta$y_centroid[i]
  annot.df <- AIH.pat3.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat3.distances <- data.frame(cell_id=AIH.pat3.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat3.distances) <- paste0("AIH.pat3_",AIH.pat3.distances$cell_id)

## AIH.pat4

AIH.pat4.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat4",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat4.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat4",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat4.meta)){
  x.pos.cell <- AIH.pat4.meta$x_centroid[i]
  y.pos.cell <- AIH.pat4.meta$y_centroid[i]
  annot.df <- AIH.pat4.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat4.distances <- data.frame(cell_id=AIH.pat4.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat4.distances) <- paste0("AIH.pat4_",AIH.pat4.distances$cell_id)

## AIH.pat5

AIH.pat5.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat5",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat5.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat5",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat5.meta)){
  x.pos.cell <- AIH.pat5.meta$x_centroid[i]
  y.pos.cell <- AIH.pat5.meta$y_centroid[i]
  annot.df <- AIH.pat5.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat5.distances <- data.frame(cell_id=AIH.pat5.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat5.distances) <- paste0("AIH.pat5_",AIH.pat5.distances$cell_id)

## AIH.pat6

AIH.pat6.TNF.pos.cells.df <- AIH.combined.integrated@meta.data[AIH.combined.integrated$TNF.status=="TNF.pos"&AIH.combined.integrated$sample=="AIH.pat6",colnames(AIH.combined.integrated@meta.data) %in% c("cell_id","x_centroid","y_centroid")]
AIH.pat6.meta <- AIH.combined.integrated@meta.data[AIH.combined.integrated$sample=="AIH.pat6",]
dist.to.closest.TNF <- numeric()
for (i in 1:nrow(AIH.pat6.meta)){
  x.pos.cell <- AIH.pat6.meta$x_centroid[i]
  y.pos.cell <- AIH.pat6.meta$y_centroid[i]
  annot.df <- AIH.pat6.TNF.pos.cells.df %>% group_by(cell_id) %>% summarise(dist.to.cell = sqrt((x_centroid - x.pos.cell)^2 + (y_centroid - y.pos.cell)^2))
  dist.to.closest.TNF[i] <- min(annot.df$dist.to.cell)
}
AIH.pat6.distances <- data.frame(cell_id=AIH.pat6.meta$cell_id,dist.to.closest.TNF=dist.to.closest.TNF)
rownames(AIH.pat6.distances) <- paste0("AIH.pat6_",AIH.pat6.distances$cell_id)

## 2. Add distances to integrated object

AIH.distances.all <- do.call(rbind,list(AIH.pat1.distances,AIH.pat2.1.distances,AIH.pat2.2.distances,AIH.pat3.distances,AIH.pat4.distances,AIH.pat5.distances,AIH.pat6.distances))
AIH.combined.integrated <- AddMetaData(AIH.combined.integrated,AIH.distances.all)
AIH.combined.integrated$TNFRSF1A.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="TNFRSF1A",]
AIH.combined.integrated$TNFRSF1B.expr <- AIH.combined.integrated[["Xenium"]]$counts[rownames(AIH.combined.integrated[["Xenium"]]$counts)=="TNFRSF1B",]


AIH.hepatocytes <- subset(AIH.combined.integrated, subset = celltype %in% c("Hepatocytes.central","Hepatocytes.portal","Hepatocytes.mixed","Hepatocytes.inflammatory"))
AIH.hepatocytes <- PrepSCTFindMarkers(AIH.hepatocytes)
dist.cutoff <- 20
AIH.hepatocytes$d20.categories <- ifelse(AIH.hepatocytes$dist.to.closest.TNF<=20,ifelse(AIH.hepatocytes$TNFRSF1A.expr>0|AIH.hepatocytes$TNFRSF1B.expr>0,"close.and.TNFRpos","close.and.TNFRneg"),ifelse(AIH.hepatocytes$TNFRSF1A.expr>0|AIH.hepatocytes$TNFRSF1B.expr>0,"far.and.TNFRpos","far.and.TNFRneg"))
AIH.markers.d20.close <- FindMarkers(AIH.hepatocytes,group.by = "d20.categories",ident.1 = "close.and.TNFRpos",ident.2 = "close.and.TNFRneg",logfc.threshold = 0,min.pct = 0,min.cells.feature = 0,min.cells.group = 0)


## DEG d20 TNFR+ vs TNFR- only Adh genes

hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))

hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )

hsa_pathways <- keggList("pathway", "hsa") %>% 
  tibble(pathway = names(.), description = .)

hsa_pathways$parent <- hsa_pathways$grand <- NA
for(i in 1:nrow(hsa_pathways)){
  #message(path)
  tryCatch({
    KeggMeta <- keggGet(hsa_pathways$pathway[i])[[1]]$CLASS
    if(is.null(KeggMeta)) next
  },
  error = function(e){
    message(e)
    next
  })
  parents = as.vector(unlist(strsplit(KeggMeta,";")))
  hsa_pathways$parent[i] <-  parents[2]
  hsa_pathways$grand[i] <- parents[1]
}

hsa_kegg_anno$pathway_ <- str_split_i(hsa_kegg_anno$pathway,pattern = ":", i=2)
hsa_kegg_anno$parent <- hsa_pathways$parent[match(hsa_kegg_anno$pathway_, hsa_pathways$pathway)]
hsa_kegg_anno$grand <- hsa_pathways$grand[match(hsa_kegg_anno$pathway_, hsa_pathways$pathway)]
hsa_kegg_anno$description <- hsa_pathways$description[match(hsa_kegg_anno$pathway_, hsa_pathways$pathway)]

AdhGenes <- hsa_kegg_anno %>% 
  subset(pathway_ == "hsa04514") %>% dplyr::select(columns = "symbol") %>% unlist %>% unname

AIH.markers.d20.close <- FindMarkers(AIH.hepatocytes,group.by = "d20.categories",ident.1 = "close.and.TNFRpos",ident.2 = "close.and.TNFRneg",logfc.threshold = 0,min.pct = 0,min.cells.feature = 0,min.cells.group = 0)
AIH.markers.d20.close$gene <- rownames(AIH.markers.d20.close)

AIH.markers.d20.close <- AIH.markers.d20.close[AIH.markers.d20.close$gene %in% AdhGenes,]

pdf(paste0(working.dir,"/results/AIH_hepat_20nm_TNFRpos_vs_TNFRneg.pdf"),width = 5, height=5)
ggplot(AIH.markers.d20.close, aes(x=avg_log2FC, y=-log10(p_val_adj), color=ifelse(p_val_adj < 0.05, "sig", "not"))) + 
  geom_point() +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 10^-0.05, linetype="dotted") +
  geom_text_repel(data= subset(AIH.markers.d20.close, p_val_adj < 0.05), aes(label=gene), color="black",point.padding = 0.25) + 
  scale_color_manual(values = c("lightsteelblue4","#F39F07")) + 
  theme_paper() + theme(legend.position = "none")
dev.off()