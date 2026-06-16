## Figure 2 ##

library(Seurat)
library(dplyr)
library(ggplot2)

CD4.combined <- readRDS(file='./Integration/RDS/12_SCS_SNS_CD4_figures.rds')
CD8.combined <- readRDS(file='./Integration/RDS/13_SCS_SNS_CD8_figures.rds')
Tcells.all.combined.subset <- readRDS(file = "./Integration/RDS/9_PBMC_liver_Tcell_clustered.rds")



# Fig 3A -----------------------------------------------------------------------

DefaultAssay(CD4.combined) <- "RNA"
cytokines.heatmap <- c("TNF","TNFSF10","FASLG",
                       "TBX21","IFNG","GATA3",
                       "IL4","IL13",
                       "RORC","IL17A","IL17F",
                       "BCL6","IL21",
                       "FOXP3","MAF","PRDM1","IL10",
                       "GZMA","GZMB","GZMH","GZMK","GZMM","PRF1")
CD4.combined <- ScaleData(CD4.combined,features = cytokines.heatmap)

avgexp.cytokine <- as.data.frame(AverageExpression(CD4.combined,features = cytokines.heatmap,assays = "RNA",slot = "scale.data")$RNA)
avgexp.cytokine.heatmap <- avgexp.cytokine[rownames(avgexp.cytokine) %in% cytokines.heatmap,] %>% rownames_to_column(var = "cytokine") %>% gather(-cytokine,key = "cell.type",value = "avg.expression")
avgexp.cytokine.heatmap$cell.type = factor(avgexp.cytokine.heatmap$cell.type, levels = c("Naive", "TCM", "TEM-a", "TEM-b", "TRM1", "TFH", "TR1", "TREG", "TREG-naive"))
avgexp.cytokine.heatmap$cytokine = factor(avgexp.cytokine.heatmap$cytokine, levels = rev(cytokines.heatmap))

pdf("SCS_SNS_CD4_heatmap_cytokines.pdf", width = 10, height = 10)
ggplot(avgexp.cytokine.heatmap, aes(x = cell.type, y = cytokine, fill = avg.expression)) +
  geom_tile(color = "white") +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradientn(colors = c("#f7fbff", "#deebf7", "#4292c6", "#2171b5", "#08519c", "#08306b"))  
dev.off()

# Fig 3B -----------------------------------------------------------------------

CD4.combined$CTaa.patient <- paste(CD4.combined$CTaa, CD4.combined$orig.ident, sep=".")
head(CD4.combined)

clonotype.df <- data.frame(orig.ident = CD4.combined$orig.ident, protocol = CD4.combined$protocol, CTaa.patient = CD4.combined$CTaa.patient)

all.clonotypes <- unique(clonotype.df$CTaa.patient)
by.clonotype <- data.frame(orig.ident=character(0),CTaa.patient=character(0),CTaa.freq=numeric(0))
for (i in 1:length(all.clonotypes)){
  clonotype.df.i <- clonotype.df[clonotype.df$CTaa.patient==all.clonotypes[i],]
  clonotype.to.add <- data.frame(orig.ident=unique(clonotype.df.i$orig.ident),
                                 CTaa.patient=unique(clonotype.df.i$CTaa.patient),
                                 CTaa.freq=sum(clonotype.df.i$protocol=="SCS"))
  head(clonotype.to.add)
  by.clonotype <- rbind(by.clonotype,clonotype.to.add)
}

CD4.combined$CTaa.freq <- NA
CD4.combined$CTaa.freq <- by.clonotype$CTaa.freq[match(CD4.combined$CTaa.patient, by.clonotype$CTaa.patient)]

CD4.combined$clonotype.add <- ifelse(CD4.combined$CTaa.freq <= 1, "Single (X = 1)", 
                            ifelse (CD4.combined$CTaa.freq <= 3, "Small (1< X <= 3)", 
                                    ifelse (CD4.combined$CTaa.freq <= 5, "Medium (3< X <= 5)",
                                            "Large (5< X <= 10)")))

Idents(CD4.combined) <- "clonotype.add"

expansion_high <- WhichCells(CD4.combined, idents = c("Large (5< X <= 10)"))
expansion_medium <- WhichCells(CD4.combined, idents = c( "Medium (3< X <= 5)"))
expansion_small <- WhichCells(CD4.combined, idents = c("Small (1< X <= 3)"))

ggplot_tcr <- DimPlot(CD4.combined, group.by = "cloneType", pt.size = 1.5,
                      cells.highlight= list(expansion_high, expansion_medium, expansion_small),
                      cols.highlight = c("#fa9fb5",  "#1f91c0", "#422977"),
                      sizes.highlight = 3)
ggsave("CD4_clonal_expansion.jpg", plot = ggplot_tcr, width = 10, height = 10, dpi = 300)

# Fig 3C -----------------------------------------------------------------------

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

Liver_expand_high <- WhichCells(CD4.combined, idents = c("Liver_expand_high"))
Dual_expand_high <- WhichCells(CD4.combined, idents = c("Dual_expand_high"))

ggplot_overlap_liver_expand_high <- DimPlot(CD4.combined, pt.size = 1.5,
                                            cells.highlight= list(Liver_expand_high),
                                            cols.highlight = c("#0570b0"),
                                            sizes.highlight = 3)
ggsave("CD4_clonaloverlap_liver_expand_high.jpg", plot = ggplot_overlap_liver_expand_high, width = 10, height = 10, dpi = 300)

ggplot_overlap_dual_expand_high <- DimPlot(CD4.combined, pt.size = 1.5,
                                           cells.highlight= list(Dual_expand_high),
                                           cols.highlight = c("#810f7c"),
                                           sizes.highlight = 3)
ggsave("CD4_clonaloverlap_dual_expand_high.jpg", plot = ggplot_overlap_dual_expand_high, width = 10, height = 10, dpi = 300)

# Fig 3D -----------------------------------------------------------------------

TRM_score <- list (c("CD69",
                     "CA10",
                     "IL17F",
                     "IL2",
                     "CDHR1",
                     "IL21",
                     "IL10",
                     "IL23R",
                     "CXCL13",
                     "CXCR6",
                     "KCNK5",
                     "ITGA1",
                     "JAG2",
                     "SRGAP3",
                     "TOX2",
                     "CH25H",
                     "NEK10",
                     "TMEM200A",
                     "MYO1B",
                     "PLXDC1",
                     "IKZF3",
                     "GFOD1",
                     "CRTAM",
                     "DUSP6",
                     "RGS1",
                     "TP53I11",
                     "GFI1",
                     "IFNG",
                     "SLC7A5",
                     "GCNT4"))
CD4.combined <- AddModuleScore(CD4.combined,
                            features = TRM_score,
                            name="TRM_score")

Migratory_score <- list(c('RIPOR2',
                            'STK38',
                            'GRASP',
                            'KLF3',
                            'SAMD3',
                            'GABBR1',
                            'FRY',
                            'ARHGEF11',
                            'VIPR1',
                            'BAIAP3',
                            'MFGE8',
                            'SBK1',
                            'HAPLN3',
                            'TTC16',
                            'CX3CR1',
                            'USP46',
                            'PLXNA4',
                            'NSG1',
                            'DSEL',
                            'CNTNAP1',
                            'VSIG1',
                            'RGMB',
                            'TTYH2',
                            'EPHA4',
                            'TNFRSF11A',
                            'MUC1',
                            'CR1',
                            'E2F2',
                            'KLF2',
                            'EDA',
                            'KRT73',
                            'ZNF462',
                            'RAP1GAP2',
                            'S1PR1',
                            'NPDC1',
                            'KLF3-AS1',
                            'ISM1',
                            'TSPAN18',
                            'KCTD15',
                            'KRT72',
                            'SEMA5A',
                            'WNT7A',
                            'SOX13',
                            'FUT7',
                            'PTGDS',
                            'PI16',
                            'SEMA3G',
                            'SYT4'))

CD4.combined <- AddModuleScore(CD4.combined,
                            features = Migratory_score,
                            name="Migratory_score")
Idents(CD4.combined) <- "annotation"
levels(CD4.combined) <- c("Naive", "TCM", "TEM-a", "TEM-b", "TRM1", "TFH", "TR1", "TREG", "TREG-naive")

pdf("SCS_SNS_CD4_migration_score_vln_corrected.pdf", width = 10, height = 10)
VlnPlot(object = CD4.combined, features = "TRM_score1",
        cols = c('TCM' = '#d99f3e',
                 'TRM1' = '#003c96',
                 'TEM-b' = '#ae98b6',
                 'TREG' = '#807dba',
                 'Naive' = 'lightgrey',
                 'TR1' = '#426e66',
                 'TFH' = '#41ab5d',
                 'TEM-a' = '#c18f6e',
                 'TREG-naive' = '#dadaeb'),
        pt.size=0.1) + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", size = 1)
VlnPlot(object = CD4.combined, features = "Migratory_score1",
        cols = c('TCM' = '#d99f3e',
                 'TRM1' = '#003c96',
                 'TEM-b' = '#ae98b6',
                 'TREG' = '#807dba',
                 'Naive' = 'lightgrey',
                 'TR1' = '#426e66',
                 'TFH' = '#41ab5d',
                 'TEM-a' = '#c18f6e',
                 'TREG-naive' = '#dadaeb'),
        pt.size=0.1) + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", size = 1)
dev.off()


# Fig 3G -----------------------------------------------------------------------

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
pdf("SCS_SNS_CD8_autoaggressive_score_CD8.pdf", width = 10, height = 10)
VlnPlot(object = CD8.combined, features = "Autoaggressive_score1",
        cols = c('TEM-a' = '#bf8fc9',
                 
                 'TEM-b' = '#9ecae1',
                 'TEM-c' = '#807dba',
                 'Naive' = '#5a1209',
                 'TEM-d' = '#ed8f87',
                 'TCM-b' = '#8c564b',
                 
                 'TEM-e' = '#f5af2d',
                 'MT-high' = 'lightgrey',
                 'TEM-f' = '#bcbd22',
                 'CD4' = '#a50f15',
                 'TCM-a' = '#225ea8',
                 'TEM-g' = '#f3bd91'),
        pt.size=0.1) + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "#fb6a4a", size = 1)
dev.off()

# Fig 3H -----------------------------------------------------------------------

CD8.combined$CTaa.patient <- paste(CD8.combined$CTaa, CD8.combined$orig.ident, sep=".")
clonotype.df <- data.frame(orig.ident = CD8.combined$orig.ident, protocol = CD8.combined$protocol, CTaa.patient = CD8.combined$CTaa.patient)
all.clonotypes <- unique(clonotype.df$CTaa.patient)

by.clonotype <- data.frame(orig.ident=character(0),CTaa.patient=character(0),CTaa.freq=numeric(0))
for (i in 1:length(all.clonotypes)){
  clonotype.df.i <- clonotype.df[clonotype.df$CTaa.patient==all.clonotypes[i],]
  clonotype.to.add <- data.frame(orig.ident=unique(clonotype.df.i$orig.ident),
                                 CTaa.patient=unique(clonotype.df.i$CTaa.patient),
                                 CTaa.freq=sum(clonotype.df.i$protocol=="SCS"))
  head(clonotype.to.add)
  by.clonotype <- rbind(by.clonotype,clonotype.to.add)
}
CD8.combined$CTaa.freq <- NA
CD8.combined$CTaa.freq <- by.clonotype$CTaa.freq[match(CD8.combined$CTaa.patient, by.clonotype$CTaa.patient)]

CD8.combined$clonotype.add <- ifelse(CD8.combined$CTaa.freq <= 1, "Single (X = 1)", 
                            ifelse (CD8.combined$CTaa.freq <= 10, "Small (1< X <= 10)", 
                                    ifelse (CD8.combined$CTaa.freq <= 30, "Medium (10< X <= 30)",
                                            "Large (30< X <= 200)")))

Idents(CD8.combined) <- "clonotype.add"

expansion_high <- WhichCells(CD8.combined, idents = c("Large (30< X <= 200)"))
expansion_medium <- WhichCells(CD8.combined, idents = c( "Medium (10< X <= 30)"))
expansion_small <- WhichCells(CD8.combined, idents = c("Small (1< X <= 10)"))

ggplot_tcr <- DimPlot(CD8.combined, pt.size = 1.5,
                      cells.highlight= list(expansion_high, expansion_medium, expansion_small),
                      cols.highlight = c("#fa9fb5",  "#1f91c0", "#422977"),
                      sizes.highlight = 2)
ggsave("CD8_clonal_expansion.jpg", plot = ggplot_tcr, width = 10, height = 10, dpi = 300)

# Fig 3I -----------------------------------------------------------------------

## same plot as in Fig 3C, but capped at 70 to visualize highly expanded CD8 clones

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

Liver_expand_high <- WhichCells(CD8.combined, idents = c("Liver_expand_high"))
Dual_expand_high <- WhichCells(CD8.combined, idents = c("Dual_expand_high"))

ggplot_overlap_liver_expand_high <- DimPlot(CD8.combined, pt.size = 1.5,
                                            cells.highlight= list(Liver_expand_high),
                                            cols.highlight = c("#0570b0"),
                                            sizes.highlight = 3)
ggsave("CD8_clonaloverlap_liver_expand_high.jpg", plot = ggplot_overlap_liver_expand_high, width = 10, height = 10, dpi = 300)

ggplot_overlap_dual_expand_high <- DimPlot(CD8.combined, pt.size = 1.5,
                                           cells.highlight= list(Dual_expand_high),
                                           cols.highlight = c("#810f7c"),
                                           sizes.highlight = 3)
ggsave("CD8_clonaloverlap_dual_expand_high.jpg", plot = ggplot_overlap_dual_expand_high, width = 10, height = 10, dpi = 300)

# Fig 3J -----------------------------------------------------------------------

DefaultAssay (CD8.combined) <- "RNA"

TRM_score <- list (c("CD69",
                     "CA10",
                     "IL17A",
                     "CXCL13",
                     "SCUBE1",
                     "HASPIN",
                     "ITGA1",
                     "CXCR6",
                     "ATP8B4",
                     "CSF1",
                     "ITGAE",
                     "CPNE7",
                     "IL10",
                     "SPRY1",
                     "MCAM",
                     "RGS1",
                     "KCNQ3",
                     "DAB2IP",
                     "TRPM2",
                     "KCNK5",
                     "IL23R",
                     "PELO",
                     "COL5A1",
                     "IRF4",
                     "FSD1",
                     "IL17RE",
                     "ADAM12",
                     "CRTAM",
                     "ARHGAP18",
                     "CCR1",
                     "JAML",
                     "ICOS",
                     "TMIGD2",
                     "TP53INP1",
                     "BMF",
                     "CD9",
                     "RIMS3",
                     "DUSP6",
                     "CCR6",
                     "GZMB",
                     "ZNF683"))
CD8.combined <- AddModuleScore(CD8.combined,
                            features = TRM_score,
                            name="TRM_score")

Circulatory_score <- list(c('PXN',
                            'FLNA',
                            'CYB561',
                            'CD300A',
                            'TSPAN32',
                            'RASA3',
                            'ADGRG5',
                            'TGFBR3',
                            'SAMD3',
                            'PELI2',
                            'C11orf21',
                            'RASGRP2',
                            'SYNE1',
                            'GK5',
                            'SSX2IP',
                            'STK38',
                            'FGR',
                            'SSBP3',
                            'CFH',
                            'ADAMTS10',
                            'MTSS1',
                            'KLF3',
                            'KLF2',
                            'SVIL',
                            'CACNA2D2',
                            'RIPOR2',
                            'SBK1',
                            'PATL2',
                            'TMCC3',
                            'KIR2DS4',
                            'HPCAL4',
                            'VCL',
                            'TTC16',
                            'PDZD4',
                            'DCHS1',
                            'EBF4',
                            'OSBPL5',
                            'FZD4',
                            'GNLY',
                            'NHSL2',
                            'TSPAN18',
                            'ME3',
                            'MSX2P1',
                            'ZNF711',
                            'NSG1',
                            'FCGR3A',
                            'GPA33',
                            'COL6A2',
                            'CXCR2',
                            'TTYH2',
                            'AGPAT4',
                            'TKTL1',
                            'SELP',
                            'LILRB1',
                            'ITGAM',
                            'LOXL4',
                            'KLF3-AS1',
                            'TFCP2L1',
                            'C1orf21',
                            'SLCO4C1',
                            'NUAK1',
                            'PALLD',
                            'DNAI2',
                            'SOX13',
                            'S1PR1',
                            'SELL',
                            'PLEKHG3',
                            'ADGRG1',
                            'SPTB',
                            'ZNF365',
                            'PCDH1',
                            'NPDC1',
                            'KRT73',
                            'KRT72',
                            'ASCL2',
                            'TAFA1',
                            'SGCD',
                            'LAIR2',
                            'EFHC2',
                            'RAP1GAP2',
                            'NME8',
                            'PODN',
                            'SH3RF2',
                            'KIF19',
                            'PTGDS',
                            'EPHX4',
                            'PRSS23',
                            'KIR3DX1',
                            'CX3CR1',
                            'SLC1A7',
                            'FGFBP2',
                            'LRFN2',
                            'DGKK'))
CD8.combined <- AddModuleScore(CD8.combined,
                            features = Circulatory_score,
                            name="Migratory_score")

pdf("SCS_SNS_CD8_migration_score_vln_corrected.pdf", width = 10, height = 10)
VlnPlot(object = CD8.combined, features = "TRM_score1",
        cols = c('TEM-a' = '#bf8fc9',
                 
                 'TEM-b' = '#9ecae1',
                 'TEM-c' = '#807dba',
                 'Naive' = '#5a1209',
                 'TEM-d' = '#ed8f87',
                 'TCM-b' = '#8c564b',
                 
                 'TEM-e' = '#f5af2d',
                 'MT-high' = 'lightgrey',
                 'TEM-f' = '#bcbd22',
                 'CD4' = '#a50f15',
                 'TCM-a' = '#225ea8',
                 'TEM-g' = '#f3bd91'),
        pt.size=0.1) + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "#fb6a4a", size = 1)
VlnPlot(object = CD8.combined, features = "Migratory_score1",
        cols = c('TEM-a' = '#bf8fc9',
                 
                 'TEM-b' = '#9ecae1',
                 'TEM-c' = '#807dba',
                 'Naive' = '#5a1209',
                 'TEM-d' = '#ed8f87',
                 'TCM-b' = '#8c564b',
                 
                 'TEM-e' = '#f5af2d',
                 'MT-high' = 'lightgrey',
                 'TEM-f' = '#bcbd22',
                 'CD4' = '#a50f15',
                 'TCM-a' = '#225ea8',
                 'TEM-g' = '#f3bd91'),
        pt.size=0.1) + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "#fb6a4a", size = 1)
dev.off()

