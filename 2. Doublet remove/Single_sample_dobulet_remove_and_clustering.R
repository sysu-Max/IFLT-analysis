library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(parallel)
library(stringi)
library(data.table)
library(RColorBrewer)
library(DoubletFinder)
library(future)
library(ggsci)
library(SoupX)
library(harmony)
library(patchwork)

TOP_N <- function(x, n, pct.1 = 0.1, first = "avg_logFC", second = "p_val_adj", sig.padj = NULL, fc.threshold = c(-0.25, 0.25)){
  if(table(names(x) %in% c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))["TRUE"] == 7){
    reordered_daf <- data.frame()
    if (!is.null(pct.1) & (pct.1 > 0 & pct.1 < 1)) {
      x <- x[x$pct.1 >= pct.1, ]
    }
    if (!is.null(sig.padj) & is.numeric(sig.padj)) {
      x <- x[x$p_val_adj <= sig.padj, ]
    }
    if (length(fc.threshold) == 2) { 
      fc.threshold <- sort(fc.threshold)
      x <- x[(x$avg_logFC < fc.threshold[1] | x$avg_logFC > fc.threshold[2]), ]
    } else if (length(fc.threshold) == 1) {
      x <- x[x$avg_logFC >= fc.threshold, ]
    }
    
    for (i in unique(x$cluster)) {
      if(n > dim(x[x$cluster == i, ])[1]){
        message("FBI warning, n < ", n)
        tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first), get(second)))[1:(dim(x[x$cluster == i, ])[1]), ]
        reordered_daf <- rbind(reordered_daf, tmp)
      } else {
        tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first), get(second)))[1:n, ]
        reordered_daf <- rbind(reordered_daf, tmp)
      }
    }
  }
  else {
    message("be careful !")
  }
  return(reordered_daf)
}

##------optional------parallel environment seting--cores = 10, memory = 100G----##
options(future.globals.maxSize = 10*1000 * 1024^2)
plan("multiprocess", workers = 10)

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]

###-----------------------1. the first clustering-----------------------####
ags <- commandArgs(trailingOnly = T) # 导入lsf 脚本里面的参数
setwd(paste0("/data1/home/wanglh/CLT_IFLT/2.Analysis"))


##------------------------1). Read in all input expression matrices
# Create and setup Seurat objects for each dataset 

TenXdat <- Read10X(data.dir = paste0("/data1/home/wanglh/CLT_IFLT/3.Script", ags[1], "/outs/filtered_feature_bc_matrix/"))
TenXdat <- CreateSeuratObject(counts = TenXdat, min.cells = 0, min.features = 500, project = ags[1])

##------------------------2).remove the genes in minimal cells (>=0.1% total cells)
# if(dim(TenXdat2@assays$RNA@data)[2] <= 3000){
#   TenXdat <- TenXdat2
# } else {
#   TenXdat <- CreateSeuratObject(counts = TenXdat1, min.cells = round(dim(TenXdat2@assays$RNA@data)[2]/1000),
#                                 min.features = 500, project = ags[1])
# }

# TenXdat <- get(load(ags[1]))

##------------------------3).cells quality filtering
mito.features <- grep(pattern = "^MT-", x = rownames(x = TenXdat), value = TRUE) # 检索所有线粒体基因
percent.mito <- Matrix::colSums(x = GetAssayData(object = TenXdat, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = TenXdat, slot = 'counts')) # 分别计算每个细胞 所有线粒体基因拷贝数之和/该细胞检测到的所有拷贝数

TenXdat[["percent.mito"]] <- percent.mito # 直接引用的meta.data 矩阵里面，barcode顺序 一致

TenXdat <- subset(x = TenXdat, subset = nFeature_RNA >= 500 & nFeature_RNA < 25000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 500000) # 基因数 500-25000，每个基因拷贝数1000-500000，线粒体基因 占比<0.25
TenXdat <- NormalizeData(object = TenXdat, normalization.method = "LogNormalize", scale.factor = 1e4) # log1p(value/colSums[cell-idx] *scale_factor); log1p <- log(1+x); 反函数 expm1() <- exp()-1; log(),即ln(); 删细胞不影响normalize
# 消除每个细胞之间测序深度不一，先除以该细胞测得的n_counts再×10^4。使得 每个细胞文库大小为10000,在log10。
TenXdat <- FindVariableFeatures(object = TenXdat, nfeatures = 2000) # 寻找2000个高变基因
TenXdat <- ScaleData(object = TenXdat, features = rownames(x = TenXdat), vars.to.regress = c("nCount_RNA", "percent.mito")) # 将归一化的基因表达转换为Z Score,一个数与平均数的差再除以标准差, 同时对表达量和线粒体基因进行回归, 并创建新的矩阵scaled.data；删细胞会影响表达矩阵

TenXdat <- RunPCA(object = TenXdat, features = VariableFeatures(object = TenXdat), verbose = FALSE) # 根据高变基因PCA降维聚类， 会出现每个基因在 PC_1、PC_2等 映射的正负值，排序 最大和最小的前top基因，默认要求数据为正态分布

# RunPCA(mast_alld, features = setdiff(mast_alld@assays$RNA@var.features, row.names(mast_alld) %>% grep(pattern = "^MT-|^RPL|^RPS|^IGK|^IGV|^IGL|^IGH", v = T)) , npcs = 50, verbose = TRUE)， 在这里可以在PCA删除目标基因，从而使其不影响UMAP


##------------------------4). dims and resulations used
subset_cells <- TenXdat

dim.use <- 30
res.use <- 1

subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use) # 使用PCA 多少个 PC_1， 一般是PCA_1：PCA_30
subset_cells <- FindClusters(object = subset_cells, resolution = res.use)

#ags <- ags %>% str_split(pattern = "_", simplify = T) %>% `[`(1)

##------------------------5). Run the UMAP
subset_cells <- RunUMAP(object = subset_cells, dims = 1:dim.use, umap.method = "uwot")

##------------------------## plot the first run tSNE
png(paste0(ags, "_", dim.use, "_", res.use, "_First_Run_tSNE.png"), height = 8, width = 10, units = 'in',res = 300)
DimPlot(object = subset_cells, reduction = 'umap', label = TRUE, cols = color_used)
dev.off()

write.table(data.frame(Tissue = ags[1], genes = dim(subset_cells@assays$RNA@data)[1], cells = dim(subset_cells@assays$RNA@data)[2]),
            file = paste0(ags[1], "_before_dobuletfinder.txt"),
            sep = "\t", row.names = F, quote = F)

##------------------------2.doublet finder processing---------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn[which(max(bcmvn$BCmetric) == bcmvn$BCmetric),2] %>% as.character %>% as.numeric

##------------------------1). Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(subset_cells@active.ident) # 即idents(subset_cells), 可以换成subset_cells$seurat_clsters， 根据我们提供的cluster，人为混合双细胞
nExp_poi <- round(0.07*length(colnames(subset_cells))) # posson 分布，查表 10000个细胞，doublet rate 为0.075, 0.075乘以总细胞
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 0.075乘以总细胞，再减去人为混合的双细胞

##------------------------2). Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
subset_cells <- doubletFinder_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE)
subset_cells <- doubletFinder_v3(subset_cells, PCs = subset_cells@commands$FindNeighbors.RNA.pca$dims, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = grep("pANN", names(subset_cells@meta.data), value = T))

##------------------------3). Plot results --------------------------------------------------------------------------------------------------------------
high_of_low <- (subset_cells@meta.data[, grep("^DF\\.classifications", names(subset_cells@meta.data), value = T)] == "Singlet")+0
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 2, "DF_hi.lo"] <- "Singlet"
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 1, "DF_hi.lo"] <- "Doublet_lo"
subset_cells@meta.data[high_of_low[, 1] + high_of_low[, 2] == 0, "DF_hi.lo"] <- "Doublet_hi"

##------------------------3.remove the doublets and recluster-----------------#######
subset_cells <- subset_cells[, subset_cells@meta.data[grepl(subset_cells$DF_hi.lo, pattern = "Singlet"), ] %>% row.names()]

##------------------------. recluster-----------------#######
subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
subset_cells <- FindVariableFeatures(object = subset_cells, selection.method = 'vst', nfeatures = 2000)
subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells), vars.to.regress = c("nCount_RNA", "percent.mito"))

subset_cells <- RunPCA(object = subset_cells, features = VariableFeatures(object = subset_cells), verbose = FALSE)

subset_cells <- FindNeighbors(object = subset_cells, dims = 1:dim.use)
subset_cells <- FindClusters(object = subset_cells, resolution = res.use)

##------------------------ Run the umap
subset_cells <- RunUMAP(object = subset_cells, dims = 1:dim.use, umap.method = "uwot")

##------------------------ plot the umap
png(paste0(ags, "_", dim.use, "_", res.use, "_doublets_removed_tSNE.png"), height = 8, width = 10, units = 'in',res = 300)
DimPlot(object = subset_cells, reduction = 'umap', label = TRUE, cols = color_used)
dev.off()

write.table(data.frame(Tissue = ags, genes = dim(subset_cells@assays$RNA@data)[1], cells = dim(subset_cells@assays$RNA@data)[2]),
            file = paste0(ags[1], "_after_dobuletfinder.txt"),
            sep = "\t", row.names = F, quote = F)

##------------------------ dfin all the markers
plan("multiprocess", workers = 1)

result <- mclapply(as.numeric(levels(subset_cells@active.ident)),
                   FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL, max.cells.per.ident = 500)},
                   mc.cores = 36)
RESULT <- result

roundN <- 1
while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
  if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
    print(recalculate_clusters)
    result1 <- mclapply(recalculate_clusters,
                        FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                        mc.cores = 35)
  }
  print(roundN + 1)
  for(i in 1:length(recalculate_clusters)){
    result[c(recalculate_clusters+1)[i]] <- result1[i]
  }
}

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result))
all_markers$cluster <- rep(levels(subset_cells@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
subset_cells.markers <- all_markers

tissue.markers <- subset_cells.markers
tissue.markers %>% TOP_N(n = 50) -> top50
write.table(top50, paste0(ags, "top50.csv"), sep = ",", row.names = T, quote = F)
write.table(tissue.markers %>% TOP_N(n = 2000), paste0(ags, ".csv"), sep = ",", row.names = T, quote = F)

##------------------------ plot the heat map

tissue.markers %>% TOP_N(n = 10) -> top10
png(paste0(ags, "_", dim.use, "_", res.use,"_doublets_removed_heatmap.png"), height = 8, width = 14, units = 'in',res = 300)
DoHeatmap(object = subset_cells, features = top10$gene, size = 2) + NoLegend() +
  theme(axis.text.x = element_text(size = 0), ##control the x label of cell barcodes
        axis.text.y = element_text(size = 0) ##control the gene label
  )
dev.off()

##------------------------ plot the vlnplot

png(paste0(ags, "_Vlnplot_of_markers.png"), height = 25, width = 15, units = "in", res = 400)
VlnPlot(subset_cells,
        features = c("CD3E", "CD3D", "FCGR3A", "MS4A1", "XBP1", "CD14", "MMP2", "ACTA2", "PECAM1", "CD8A", "CD8B", "CD4", "TRDV1", "TRDV2", "C1QC", "EPCAM", "ALB"), ## active and repression marker genes
        # features = c("MKI67", "PCNA", "TYMS"),
        pt.size = 0.01,
        log = F,
        ncol = 1,
        cols = color_used)  
dev.off()


# subset_cells <- subset(subset_cells, idents = samples)
sc <-  SoupChannel(subset_cells@assays$RNA@counts, 
                   subset_cells@assays$RNA@counts, calcSoupProfile = FALSE)

toc <-  sc$toc
soupProf <-  data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
sc <- setSoupProfile(sc, soupProf)

sc <-  setClusters(sc, subset_cells$seurat_clusters %>% as.character %>% as.numeric())

# igGenes = c("ALB", "APOA2", "APOC1", "APOC3", "RBP4", "APOA1")#Z762
# igGenes = c("ALB", "APOA2", "APOC1", "APOC3", "RBP4", "APOA1", "CLPS", "PRSS1", "PRSS2", "SPINK1", "CPA1")#Z765

# useToEst <-  estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes))
# sc <-  calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst, forceAccept = TRUE)

sc <-  autoEstCont(sc)
out <-  adjustCounts(sc)

subset_cells@assays$RNA@counts <- round(out)
subset_cells <- subset(x = subset_cells, subset = nFeature_RNA >= 500 & nFeature_RNA < 25000 & percent.mito <= 0.25 & nCount_RNA >= 1000 & nCount_RNA <= 500000)
subset_cells <- NormalizeData(object = subset_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
subset_cells <- FindVariableFeatures(object = subset_cells, nfeatures = 2000)
subset_cells <- ScaleData(object = subset_cells, features = rownames(x = subset_cells), vars.to.regress = c("nCount_RNA", "percent.mito"))


# save(list = "subset_cells", file = paste0("/data1/home/wanglh/CLT_IFLT/Contaminated_gene_removed_", ags[1], ".RData"))
subset_cells[["Sample_ID"]] <- ags[1]
assign(ags[1], subset_cells)
# 
# dat <- list()
###it is so import that "list" args was used in save function.
save(list = ags[1], file = paste0("/data1/home/wanglh/CLT_IFLT/Contaminated_gene_removed_", ags[1], ".RData"))



