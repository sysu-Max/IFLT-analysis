
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(ggthemes)
library(parallel)
library(dplyr)
library(Seurat)
library(ggdendro)
library(dendextend)
library(tidyverse)
library(reshape2)
library(GSEABase)
library(SCENIC)
library(ggsci)
library(AUCell)
library(ComplexHeatmap)
# http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
# This tutorial goes through the steps in the SCENIC workflow:
#   
# 1. Building the gene regulatory network (GRN):
#   
# 2. Identify potential targets for each TF based on co-expression.
# 3. Filtering the expression matrix and running GENIE3/GRNBoost.
# 4. Formatting the targets from GENIE3/GRNBoost into co-expression modules.
# 5. Select potential direct-binding targets (regulons) based on DNA-motif analysis (RcisTarget: TF motif analysis)

setwd('~/CLT_IFLT/2.Analysis/merged/Scenic/SNP')
org <- "hgnc" # or hgnc, or dmel 
dbDir <- "~/RAW_data/1-SingleCell/3-HCA/3-analysis/8-reanalysis/All_singlet_data/" # RcisTarget databases location
myDatasetTitle <- "IFLT" # choose a name for your analysis
data(defaultDbNames)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle, nCores = 16) 
saveRDS(scenicOptions, file = paste0("scenicOptions.Rds"))

# ---load data
# 1
load(file = '~/CLT_IFLT/2.Analysis/merged/Annotation_data.RData')
donor_tsv <- openxlsx::read.xlsx('~/CLT_IFLT/2.Analysis/merged/SNP/tsv/donor_tsv.xlsx')
receptor_tsv <- openxlsx::read.xlsx('~/CLT_IFLT/2.Analysis/merged/SNP/tsv/recipient_tsv.xlsx')

subset_cells$plot_clusters <- subset_cells$new_clusters %>% as.character()
subset_cells$time_point <- subset_cells$sample.group %>% as.character()
subset_cells$time_point[subset_cells$sample.group %in% c('CLT_EP','IFLT_EP')] <- 'EP'
subset_cells$time_point[subset_cells$sample.group %in% c('CLT_PR','IFLT_PR')] <- 'PR'
subset_cells$identity <- subset_cells$sample.group %>% as.character()
subset_cells$identity[row.names(subset_cells@meta.data) %in% donor_tsv$new_barcode] <- 'donor'
subset_cells$identity[row.names(subset_cells@meta.data) %in% receptor_tsv$new_barcode] <- 'receptor'
plot_cell <- subset(subset_cells, identity %in% c('donor','receptor'))
plot_cell$identity_sample <- paste(plot_cell$identity, plot_cell$sample.group, sep = '-')

Idents(plot_cell) <- plot_cell$identity_sample
subset_cells_selected <- subset(plot_cell, cells=WhichCells(plot_cell, downsample=500))

cellInfo <- data.frame(seuratCluster = Idents(subset_cells_selected))
singleCellMatrix <- subset_cells_selected@assays$RNA@counts %>% as.matrix


genesKept <- geneFiltering(singleCellMatrix, scenicOptions = scenicOptions,
                           minCountsPerGene = 1*.01*ncol(singleCellMatrix),
                           minSamples = ncol(singleCellMatrix)*.003)

exprMat_filtered <- singleCellMatrix[genesKept, ]
dim(exprMat_filtered)

runCorrelation(exprMat_filtered, scenicOptions)
exprMat_log <- log2(singleCellMatrix + 1) 

# Run GENIE3
exprMat_filtered_log <- log2(exprMat_filtered + 1) 
runGenie3(exprMat_filtered_log, scenicOptions)


#############-----------------------------###################
scenicOptions <- readRDS("~/CLT_IFLT/2.Analysis/merged/Scenic/SNP/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 16
scenicOptions@settings$seed <- 123


runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)


cellInfo <- cellInfo


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))



##---remove extend TFs---------------####
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType[, ]), center = T, scale=T))
row.names(regulonActivity_byCellType_Scaled) <- regulonActivity_byCellType_Scaled %>% row.names %>% gsub(., pattern="_extended*", replacement = "")


cellInfos <- str_split(colnames(regulonAUC@assays@data$AUC), "_", simplify = T)[, 1]


col <- colorRampPalette(c("#2166ac", "white", "#CC3333"))(100)


anto <- data.frame(cluster = colnames(regulonActivity_byCellType_Scaled))
row.names(anto) <- colnames(regulonActivity_byCellType_Scaled)

cluster <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12))[-8][1:(cellInfo$seuratCluster %>% unique %>% length)]
names(cluster) <- anto$cluster

ancols <- list(cluster = cluster)

saveRDS(regulonAUC,file = "regulonAUC.rds")
saveRDS(regulonActivity_byCellType_Scaled,file = "regulonActivity_byCellType_Scaled.rds")

regulonActivity_byCellType_Scaled <- t(regulonActivity_byCellType_Scaled)
png("scenic_heatmap_1.png", res = 200,units = 'in',height = 10, width = 30)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled[,], #fontsize_row=3, 
                   color = colorRampPalette(c("#2166ac", "white", "#CC3333"))(100),
                   breaks = seq(-3, 3, length.out = 100),
                   treeheight_row = 0,
                   treeheight_col = 0,
                   border_color = "white",
                   cellwidth = 8,
                   cellheight = 8,
                   fontsize_row  = 4,
                   fontsize_col = 8,
                   cluster_rows = T,
                   cluster_cols = T,
                   scale = "none",
                   width = 10,
                   height = 14)
dev.off()

write.csv(data.frame(tf_s = colnames(regulonActivity_byCellType_Scaled)),file = 'tf_s.csv',row.names = F)
pdf("scenic_heatmap_1.pdf", height = 10, width = 30)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled[,], #fontsize_row=3, 
                   color = colorRampPalette(c("#2166ac", "white", "#CC3333"))(100),
                   breaks = seq(-3, 3, length.out = 100),
                   treeheight_row = 0,
                   treeheight_col = 0,
                   border_color = "white",
                   cellwidth = 8,
                   cellheight = 8,
                   fontsize_row  = 4,
                   fontsize_col = 8,
                   cluster_rows = T,
                   cluster_cols = T,
                   scale = "none",
                   width = 10,
                   height = 14)
dev.off()

tsneAUC(scenicOptions, aucType="AUC")
regulonActivity_byCellType_Scaled %>% openxlsx::write.xlsx(file = 'regulonActivity_byCellType_Scaled.xlsx', rowNames = T,
                                                           colNames =T)

