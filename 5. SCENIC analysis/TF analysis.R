library(ggplot2)
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(ggthemes)
library(dplyr)
library(Seurat)
library(reshape2)
library(GSEABase)
library(SCENIC)
library(ggsci)
library(AUCell)
library(ComplexHeatmap)
library(tidyverse)
library(pheatmap)
library(viridis)
library(pheatmap)
library(reshape2)
library(scales)
library(rlang)
library(DoubletFinder)
library(future)
library(ggsci)
library(SoupX)
library(harmony)
library(patchwork)
library(scales)
library(tidyverse)



####--scenic find clusters---
donor_tsv <- openxlsx::read.xlsx('/data1/home/wanglh/CLT_IFLT/2.Analysis/merged/SNP/tsv/donor_tsv.xlsx')
recipient_tsv <- openxlsx::read.xlsx('/data1/home/wanglh/CLT_IFLT/2.Analysis/merged/SNP/tsv/recipient_tsv.xlsx')

load(file = '/data1/home/wanglh/CLT_IFLT/2.Analysis/merged/Annotation_data.RData')
load(file = '/data1/home/wanglh/CLT_IFLT/2.Analysis/merged/Annotation_data.RData')
subset_cells$plot_clusters <- subset_cells$new_clusters %>% as.character()
subset_cells$identity <- subset_cells$sample.group %>% as.character()
subset_cells$identity[row.names(subset_cells@meta.data) %in% donor_tsv$new_barcode] <- 'donor'
subset_cells$identity[row.names(subset_cells@meta.data) %in% recipient_tsv$new_barcode] <- 'recipient'
subset_cells <- subset(subset_cells, identity %in% c('donor','recipient')) ###删除错分的细胞

subset_cells$identity_sample <- paste(subset_cells$identity,merged_CD4_7$sample.group, sep = '-')

setwd('/data1/home/wanglh/CLT_IFLT/2.Analysis/merged/Scenic/SNP')
subset_cells_selected <- subset(subset_cells,  plot_clusters == 'target_clusters') # selected target clusters
Idents(subset_cells_selected) <- subset_cells_selected$identity_sample
scenic_cells <- subset(subset_cells_selected, cells=WhichCells(subset_cells_selected, downsample=500))

####
regulonAUC <- readRDS(file = '/data1/home/wanglh/CLT_IFLT/2.Analysis/merged/Scenic/SNP/regulonAUC.rds')
scenic_cells@assays$RNA@data <- regulonAUC@assays@data@listData$AUC
genes_modified <- scenic_cells %>% row.names %>% gsub(., pattern="_extended*", replacement = "")
rownames(scenic_cells@assays$RNA@data) <- genes_modified

Idents(scenic_cells) <- scenic_cells$identity_sample

scenic_cells@active.ident <- scenic_cells@active.ident[colnames(scenic_cells@assays$RNA@data)]
scenic_cells@meta.data <- scenic_cells@meta.data[colnames(scenic_cells@assays$RNA@data),]

Find_markers <- FindAllMarkers(scenic_cells, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.01)

top50 <- Find_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% openxlsx::write.xlsx(file = 'top50_scenic_1.xlsx',rowNames = T, colNames = T)

