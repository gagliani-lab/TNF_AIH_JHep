## Extended Figure 5 ##

library(Seurat)
library(dplyr)
library(ggplot2)
library(nlme)
library(multcomp)

CD8.combined <- readRDS(file='./Integration/RDS/13_SCS_SNS_CD8_figures.rds')
Tcells.all.combined.subset <- readRDS(file = "./Integration/RDS/9_PBMC_liver_Tcell_clustered.rds")


# ED Fig 5A --------------------------------------------------------------------

DefaultAssay (CD8.combined) <- "RNA"
feature_plots <- FeaturePlot(CD8.combined, features = c("TNF", "IFNG", "PRF1", "GZMB", "GZMH", "FASLG", "CD69", "CXCR6"),
                             cols = c("lightgrey", "#2171b5"),
                             order = T,
                             slot = "data",
                             min.cutoff = "q5", max.cutoff = "q95",
                             reduction = "umap",
                             pt.size = 1.5,
                             combine= FALSE)

if (is.list(feature_plots)) {
  # Save each plot separately
  for (i in seq_along(feature_plots)) {
    ggsave(paste0("feature_plot_CD8_RNA_", feature_plots[[i]]@labels$title, ".jpg"), plot = feature_plots[[i]], width = 10, height = 10, dpi = 300)
  }
} else {
  # If it's a single plot, save it directly
  ggsave(paste0("feature_plot_CD8_RNA_", feature_plots@labels$title, ".jpg"), plot = feature_plots, width = 10, height = 10, dpi = 300)
}

DefaultAssay (CD8.combined) <- "CITE"
feature_plots_2 <- FeaturePlot(CD8.combined, features = c("ADT-CD69"),
                               cols = c("lightgrey", "#3f007d"),
                               order = T,
                               slot = "data",
                               min.cutoff = "q5", max.cutoff = "q95",
                               reduction = "umap",
                               pt.size = 2,
                               combine= FALSE)

if (is.list(feature_plots_2)) {
  # Save each plot separately
  for (i in seq_along(feature_plots_2)) {
    ggsave(paste0("feature_plot_CD8_CITE_", feature_plots_2[[i]]@labels$title, ".jpg"), plot = feature_plots_2[[i]], width = 10, height = 10, dpi = 300)
  }
} else {
  # If it's a single plot, save it directly
  ggsave(paste0("feature_plot_CD8_CITE_", feature_plots_2@labels$title, ".jpg"), plot = feature_plots_2, width = 10, height = 10, dpi = 300)
}


# ED Fig 5B --------------------------------------------------------------------

Auto_score <- list (c("RUNX3",
                      "JUN",
                      "RGS1",
                      "DUSP2",
                      "NR4A2",
                      "PDCD1",
                      "CXCR6",
                      "TIGIT",
                      "ZBTB38",
                      "GZMK",
                      "GZMA",
                      "CD74",
                      "FYN",
                      "IFNGR1",
                      "TOX",
                      "GABARAPL1",
                      "IFNG",
                      "ALOX5AP",
                      "CCL5",
                      "IKZF3",
                      "CST7",
                      "GZMM",
                      "ZFP36",
                      "CD8A"))

DefaultAssay(CD8.combined) <- "RNA"
CD8.combined<- AddModuleScore(CD8.combined,
                              features = Auto_score,
                              name="Autoaggressive_score")

## extract data to calculate statistics
data.autoaggressive <- CD8.combined@meta.data[,colnames(CD8.combined@meta.data) %in% c("orig.ident","annotation","Autoaggressive_score1")]

## remove all cells not TCM or TEM
data.autoaggressive.filtered <- data.autoaggressive[grepl("TEM",data.autoaggressive$annotation)|grepl("TCM",data.autoaggressive$annotation),]

data.autoaggressive$annotation <- factor(data.autoaggressive$annotation,levels = c("Naive", "TCM-a", "TCM-b", "TEM-a", "TEM-b", "TEM-c", "TEM-d", "TEM-e", "TEM-f", "TEM-g", "MT-high", "CD4"))

data.autoaggressive.avg <- data.autoaggressive %>% group_by(orig.ident,annotation) %>% summarise(av.autoaggressive.score = mean(Autoaggressive_score1))
data.autoaggressive.avg$protocol <- ifelse(grepl("NUC",data.autoaggressive.avg$orig.ident),"SNS","SCS")

data.autoaggressive$annotation <- factor(data.autoaggressive$annotation,levels = c("Naive", "TCM-a", "TCM-b", "TEM-a", "TEM-b", "TEM-c", "TEM-d", "TEM-e", "TEM-f", "TEM-g", "MT-high", "CD4"))

pdf("C:/Yang_AIH_Final_Version/Integration/results/251106_CD8_by_donor_larger_CD8A.pdf",width = 7, height = 5)
p1 <- ggplot(data.autoaggressive.avg,aes(x=annotation,y=av.autoaggressive.score,fill=annotation)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values = c('TEM-b' = '#9ecae1',
                               'TEM-c' = '#807dba',
                               'Naive' = '#5a1209',
                               'TEM-d' = '#ed8f87',
                               'TCM-b' = '#8c564b',
                               'TEM-a' = '#bf8fc9',    
                               'TEM-e' = '#f5af2d',
                               'MT-high' = 'lightgrey',
                               'TEM-f' = '#bcbd22',
                               'CD4' = '#a50f15',
                               'TCM-a' = '#225ea8', 
                               'TEM-g' = '#f3bd91')) +
  geom_jitter(aes(shape = protocol),size = 2.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p1)
dev.off()

## Statistics:

## combine the two TCM clusters for statistical testing only
data.autoaggressive.filtered$clusters <- as.character(data.autoaggressive.filtered$annotation)
data.autoaggressive.filtered$clusters[grepl("TCM",data.autoaggressive.filtered$clusters)] <- "TCM"
data.autoaggressive.filtered$clusters <- factor(data.autoaggressive.filtered$clusters,levels=sort(unique(data.autoaggressive.filtered$clusters)))

data.autoaggressive.filtered.avg <- data.autoaggressive.filtered %>% group_by(orig.ident,clusters) %>% summarise(av.autoaggressive.score = mean(Autoaggressive_score1))
head(data.autoaggressive.filtered.avg)

model <- lme(av.autoaggressive.score ~ clusters, random=~1|orig.ident, data = data.autoaggressive.filtered.avg)
summary(model)
post.hoc <- glht(model, linfct = mcp(clusters = 'Dunnett'))
summary(post.hoc)


# ED Fig 5C --------------------------------------------------------------------

df_freq <- CD8.combined@meta.data

proportions.clonotypes <- as.data.frame(prop.table(table(df_freq$annotation, df_freq$clonotype.add,useNA="ifany"),margin = 1))
colnames(proportions.clonotypes) <- c("annotation","clone.type","proportion")
proportions.clonotypes$percentage <- proportions.clonotypes$proportion*100
proportions.clonotypes$clone.type <- factor(proportions.clonotypes$clone.type, levels=rev(c('Single (X = 1)','Small (1< X <= 10)','Medium (10< X <= 30)','Large (30< X <= 200)')))
proportions.clonotypes$annotation <- factor(proportions.clonotypes$annotation, levels =c("Naive", "TCM-a", "TCM-b", "TEM-a", "TEM-b", "TEM-c", "TEM-d", "TEM-e", "TEM-f", "TEM-g", "MT-high", "CD4"))

pdf("SCS_SNS_CD8_clone_types_proportion.pdf", width = 10, height = 10)
ggplot(proportions.clonotypes, aes(x=annotation, y=percentage, fill=clone.type)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = rev(c("lightgrey","#fa9fb5",  "#1f91c0", "#422977"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()



# ED Fig 5D --------------------------------------------------------------------

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


CD3.by.clonotype$n.liver.capped.2 <- ifelse(CD3.by.clonotype$n.liver>=70,70,CD3.by.clonotype$n.liver)
CD3.by.clonotype$n.blood.capped.2 <- ifelse(CD3.by.clonotype$n.blood>=70,70,CD3.by.clonotype$n.blood)

pdf("CD3_clonotypes_by_number_cap70.pdf", width = 10, height = 10)
ggplot(CD3.by.clonotype, aes(x=n.liver.capped.2,y=n.blood.capped.2)) +
  geom_point(pch = 16, position = position_jitter(),size=5,color=rgb(1,0,0,0.4)) +
  theme_classic() +
  scale_y_continuous(limits = c(-0.5,75),breaks = seq(0,70,10)) +
  scale_x_continuous(limits = c(-0.5,75),breaks = seq(0,70,10))
dev.off()

## Overlay expanded clonotypes on CD8 UMAP

clone.overlap.meta <- data.frame(barcode = rownames(CD3_df), cloneoverlap = CD3_df[["cloneoverlap"]])
rownames(cloneoverlap_table) <- cloneoverlap_table$barcode
cloneoverlap_table$barcode <- NULL

# Filter cloneoverlap_table while preserving row names
filt.clone.overlap.meta <- clone.overlap.meta[rownames(clone.overlap.meta) %in% colnames(CD8.combined), , drop = FALSE]

CD8.combined <- AddMetaData(object = CD8.combined, metadata = filt.clone.overlap.meta, col.name = colnames(filt.clone.overlap.meta))


CD8.combined@meta.data <- CD8.combined@meta.data %>%
  mutate(
    Freq_overlap = case_when(
      CTaa.freq > 30 & cloneoverlap == "Liver" ~ "Liver_expand_high",
      CTaa.freq > 30 & cloneoverlap == "dual" ~ "Dual_expand_high", 
      CTaa.freq >1 & CTaa.freq <= 30 & cloneoverlap == "Liver" ~ "Liver_expand_low",
      CTaa.freq >1 & CTaa.freq <= 30 & cloneoverlap == "dual" ~ "Dual_expand_low",      
      CTaa.freq == 1 & cloneoverlap == "Liver" ~ "Liver_1",
      CTaa.freq == 1 & cloneoverlap == "dual" ~ "Dual_1",
      TRUE ~ NA_character_
    )
  )
Idents(CD8.combined) <- "Freq_overlap"

Liver_expand_low <- WhichCells(CD8.combined, idents = c("Liver_expand_low"))
Dual_expand_low <- WhichCells(CD8.combined, idents = c("Dual_expand_low"))

ggplot_overlap_liver_expand_low <- DimPlot(CD8.combined, pt.size = 1.5,
                                           cells.highlight= list(Liver_expand_low),
                                           cols.highlight = c("#0570b0"),
                                           sizes.highlight = 3)
ggsave("CD8_clonaloverlap_liver_expand_low.jpg", plot = ggplot_overlap_liver_expand_low, width = 10, height = 10, dpi = 300)

ggplot_overlap_dual_expand_low <- DimPlot(CD8.combined, pt.size = 1.5,
                                          cells.highlight= list(Dual_expand_low, Dual_1),
                                          cols.highlight = c("#810f7c"),
                                          sizes.highlight = 3)
ggsave("CD8_clonaloverlap_dual_expand_low.jpg", plot = ggplot_overlap_dual_expand_low, width = 10, height = 10, dpi = 300)