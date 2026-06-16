## Extended Figure 4 ##

library(Seurat)
library(dplyr)
library(ggplot2)


CD4.combined <- readRDS(file='./Integration/RDS/12_SCS_SNS_CD4_figures.rds')
Tcells.all.combined.subset <- readRDS(file = "./Integration/RDS/9_PBMC_liver_Tcell_clustered.rds")


# ED Fig 4A --------------------------------------------------------------------
df_freq <- CD4.combined@meta.data

proportions.clonotypes <- as.data.frame(prop.table(table(df_freq$annotation, df_freq$clonotype.add,useNA="ifany"),margin = 1))
colnames(proportions.clonotypes) <- c("annotation","clone.type","proportion")
proportions.clonotypes$percentage <- proportions.clonotypes$proportion*100
proportions.clonotypes$clone.type <- factor(proportions.clonotypes$clone.type, levels=rev(c('Single (X = 1)','Small (1< X <= 3)','Medium (3< X <= 5)','Large (5< X <= 10)')))
proportions.clonotypes$annotation <- factor(proportions.clonotypes$annotation, levels =c("Naive", "TCM", "TEM-a", "TEM-b", "TRM1", "TFH", "TR1", "TREG", "TREG-naive"))

pdf("SCS_SNS_CD4_clone_types_proportion.pdf", width = 10, height = 10)
ggplot(proportions.clonotypes, aes(x=annotation, y=percentage, fill=clone.type)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = rev(c("lightgrey","#fa9fb5",  "#1f91c0", "#422977"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# ED Fig 4C --------------------------------------------------------------------

ggplot1 <- DimPlot(Tcells.all.combined.subset, 
                   label = FALSE, 
                   group.by = "tissue", 
                   pt.size=0.5, 
                   cols=c('#cb5148','#6a9bcb'))
ggsave("PBMC_Liver_091315_Tcell_umap.jpg", plot = ggplot1, width = 10, height = 10, dpi = 300)

TRM1.barcodes <- read.delim("./Integration/Table/TRM1_barcode.txt",header = FALSE)$V1

Tcells.all.combined.subset$TRM1.status <- ifelse(rownames(Tcells.all.combined.subset@meta.data) %in% TRM1.barcodes,"Positive","Negative")

Idents(Tcells.all.combined.subset) <- "TRM1.status"
Positive <- WhichCells(Tcells.all.combined.subset, idents = c("Positive"))

ggplot2 <- DimPlot(Tcells.all.combined.subset, pt.size = 0.5,
                   cells.highlight= list(Positive),
                   cols.highlight = c("#003c96"),
                   sizes.highlight = 0.5)
ggsave("PBMC_Liver_091315_TRM1_umap.jpg", plot = ggplot2, width = 10, height = 10, dpi = 300)

aaCD8.barcodes <- read.delim("./Integration/Table/aaCD8_barcode.txt",header = FALSE)$V1
Tcells.all.combined.subset$aaCD8.status <- ifelse(rownames(Tcells.all.combined.subset@meta.data) %in% aaCD8.barcodes,"Positive","Negative")

Idents(Tcells.all.combined.subset) <- "aaCD8.status"
Positive <- WhichCells(Tcells.all.combined.subset, idents = c("Positive"))

ggplot3 <- DimPlot(Tcells.all.combined.subset, pt.size = 0.5,
                   cells.highlight= list(Positive),
                   cols.highlight = c("#6a51a3"),
                   sizes.highlight = 0.5)
ggsave("PBMC_Liver_091315_aaCD8_umap.jpg", plot = ggplot3, width = 10, height = 10, dpi = 300)


# ED Fig 4D --------------------------------------------------------------------

CD3_df <- Tcells.all.combined.subset@meta.data

# Loop through each row
for (i in 1:nrow(CD3_df)) {
  # Check if there is a duplicate of the current row based on CTaa
  if (sum(CD3_df$CTaa == CD3_df$CTaa[i]) > 1) {
    # Check if any of the duplicates have a different tissue value
    if (any(CD3_df$tissue[CD3_df$CTaa == CD3_df$CTaa[i]] != CD3_df$tissue[i])) {
      # Assign "dual" if there is a duplicate with a different tissue value
      CD3_df$cloneoverlap[i] <- "dual"
    } else {
      # Otherwise, assign the same value as tissue
      CD3_df$cloneoverlap[i] <- CD3_df$tissue[i]
    }
  } else {
    # If there are no duplicates based on CTaa, assign the same value as tissue
    CD3_df$cloneoverlap[i] <- CD3_df$tissue[i]
  }
}

CD3_df$patient <- ifelse(grepl("A09",CD3_df$orig.ident),"A09",ifelse(grepl("A15",CD3_df$orig.ident),"A15","A13"))

CD3_df$CTaa.patient <- paste(CD3_df$CTaa, CD3_df$patient, sep=".")

clonotype.df.CD3 <- CD3_df[,colnames(CD3_df) %in% c("tissue","cloneoverlap","patient","CTaa.patient")]

all.clonotypes <- unique(clonotype.df.CD3$CTaa.patient)

CD3.by.clonotype <- data.frame(patient=character(0),
                               cloneoverlap=character(0),
                               CTaa.patient=character(0),
                               n.blood=numeric(0),
                               n.liver=numeric(0))

for (i in 1:length(all.clonotypes)){
  clonotype.df.i <- clonotype.df.CD3[clonotype.df.CD3$CTaa.patient==all.clonotypes[i],]
  clonotype.to.add <- data.frame(patient=unique(clonotype.df.i$patient),
                                 cloneoverlap=unique(clonotype.df.i$cloneoverlap),
                                 CTaa.patient=unique(clonotype.df.i$CTaa.patient),
                                 n.blood=sum(clonotype.df.i$tissue=="PBMC"),
                                 n.liver=sum(clonotype.df.i$tissue=="Liver"))
  head(clonotype.to.add)
  CD3.by.clonotype <- rbind(CD3.by.clonotype,clonotype.to.add)
}

CD3.by.clonotype$n.liver.capped <- ifelse(CD3.by.clonotype$n.liver>=20,20,CD3.by.clonotype$n.liver)
CD3.by.clonotype$n.blood.capped <- ifelse(CD3.by.clonotype$n.blood>=20,20,CD3.by.clonotype$n.blood)

pdf("CD3_clonotypes_by_number_cap20.pdf", width = 10, height = 10)
ggplot(CD3.by.clonotype, aes(x=n.liver.capped,y=n.blood.capped)) +
  geom_point(pch = 16, position = position_jitter(),size=5,color=rgb(1,0,0,0.4)) +
  theme_classic() +
  scale_y_continuous(limits = c(-0.5,20.5),breaks = seq(0,20,2)) +
  scale_x_continuous(limits = c(-0.5,20.5),breaks = seq(0,20,2))
dev.off()

## Overlay expanded clonotypes on CD4 UMAP

clone.overlap.meta <- data.frame(barcode = rownames(CD3_df), cloneoverlap = CD3_df[["cloneoverlap"]])
rownames(cloneoverlap_table) <- cloneoverlap_table$barcode
cloneoverlap_table$barcode <- NULL

# Filter cloneoverlap_table while preserving row names
filt.clone.overlap.meta <- clone.overlap.meta[rownames(clone.overlap.meta) %in% colnames(CD4.combined), , drop = FALSE]

CD4.combined <- AddMetaData(object = CD4.combined, metadata = filt.clone.overlap.meta, col.name = colnames(filt.clone.overlap.meta))

CD4.combined@meta.data <- CD4.combined@meta.data %>%
  mutate(
    Freq_overlap = case_when(
      CTaa.freq > 5 & cloneoverlap == "Liver" ~ "Liver_expand_high",
      CTaa.freq > 5 & cloneoverlap == "dual" ~ "Dual_expand_high", 
      CTaa.freq >1 & CTaa.freq <= 5 & cloneoverlap == "Liver" ~ "Liver_expand_low",
      CTaa.freq >1 & CTaa.freq <= 5 & cloneoverlap == "dual" ~ "Dual_expand_low",       
      CTaa.freq == 1 & cloneoverlap == "Liver" ~ "Liver_1",
      CTaa.freq == 1 & cloneoverlap == "dual" ~ "Dual_1",
      TRUE ~ NA_character_
    )
  )

Idents(CD4.combined) <- "Freq_overlap"

Liver_expand_low <- WhichCells(CD4.combined, idents = c("Liver_expand_low"))
Dual_expand_low <- WhichCells(CD4.combined, idents = c("Dual_expand_low"))

ggplot_overlap_liver_expand_low <- DimPlot(CD4.combined, pt.size = 1.5,
                                            cells.highlight= list(Liver_expand_low),
                                            cols.highlight = c("#0570b0"),
                                            sizes.highlight = 3)
ggsave("CD4_clonaloverlap_liver_expand_low.jpg", plot = ggplot_overlap_liver_expand_low, width = 10, height = 10, dpi = 300)

ggplot_overlap_dual_expand_low <- DimPlot(CD4.combined, pt.size = 1.5,
                                           cells.highlight= list(Dual_expand_low),
                                           cols.highlight = c("#810f7c"),
                                           sizes.highlight = 3)
ggsave("CD4_clonaloverlap_dual_expand_low.jpg", plot = ggplot_overlap_dual_expand_low, width = 10, height = 10, dpi = 300)

# ED Fig 4E --------------------------------------------------------------------

DefaultAssay (CD4.combined) <- "RNA"
feature_plots <- FeaturePlot(CD4.combined, features = c("CD69", "CXCR6", "ITGA1", "TIGIT"),
                             cols = c("lightgrey", "#2171b5"),
                             order = T,
                             slot = "data",
                             min.cutoff = "q5", max.cutoff = "q95",
                             reduction = "umap",
                             pt.size = 2,
                             combine= FALSE)
feature_plots

if (is.list(feature_plots)) {
  # Save each plot separately
  for (i in seq_along(feature_plots)) {
    ggsave(paste0("feature_plot_CD4_RNA_", feature_plots[[i]]@labels$title, ".jpg"), plot = feature_plots[[i]], width = 10, height = 10, dpi = 300)
  }
} else {
  # If it's a single plot, save it directly
  ggsave(paste0("feature_plot_CD4_RNA_", feature_plots@labels$title, ".jpg"), plot = feature_plots, width = 10, height = 10, dpi = 300)
}

DefaultAssay (CD4.combined) <- "CITE"
feature_plots_2 <- FeaturePlot(CD4.combined, features = c("ADT-CD69"),
                               cols = c("lightgrey", "#3f007d"),
                               order = T,
                               slot = "data",
                               min.cutoff = "q5", max.cutoff = "q95",
                               reduction = "umap",
                               pt.size = 2,
                               combine= FALSE)
feature_plots_2

if (is.list(feature_plots_2)) {
  # Save each plot separately
  for (i in seq_along(feature_plots_2)) {
    ggsave(paste0("feature_plot_CD4_CITE_", feature_plots_2[[i]]@labels$title, ".jpg"), plot = feature_plots_2[[i]], width = 10, height = 10, dpi = 300)
  }
} else {
  # If it's a single plot, save it directly
  ggsave(paste0("feature_plot_CD4_CITE_", feature_plots_2@labels$title, ".jpg"), plot = feature_plots_2, width = 10, height = 10, dpi = 300)
}