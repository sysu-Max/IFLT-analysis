library(scRepertoire)
library(Seurat)
library(limma)
library(tidyverse)
library(scales)
library(ggsci)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(RColorBrewer)
library("viridis")  
library(parallel)
library(pbapply)
library(dplyr)
library(migest)
library(circlize)




  
####----load data-----
###TCR
D1_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D1/filtered_contig_annotations.csv')
D2_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D2/filtered_contig_annotations.csv')
D3_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D3/filtered_contig_annotations.csv')
D4_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D4/filtered_contig_annotations.csv')
D5_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D5/filtered_contig_annotations.csv')
D6_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D6/filtered_contig_annotations.csv')
D7_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D7/filtered_contig_annotations.csv')
D8_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D8/filtered_contig_annotations.csv')
D9_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D9/filtered_contig_annotations.csv')
D10_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D10/filtered_contig_annotations.csv')
D11_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D11/filtered_contig_annotations.csv')
D12_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D12/filtered_contig_annotations.csv')
D15_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D15/filtered_contig_annotations.csv')
D16_TCR <- read.csv('D:/IFLT-CLT/2.Analysis/integrated/TCR/D16/filtered_contig_annotations.csv')



###--barcode celltype--
D1_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD1.csv',header = T)
D1_barcode$cell_barcode <- str_split_fixed(D1_barcode$cell_barcode,pattern = '_',n = 2)[,2]
D2_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD2.csv',header = T)
D2_barcode$cell_barcode <- str_split_fixed(D2_barcode$cell_barcode,pattern = '_',n = 2)[,2]

D3_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD3.csv',header = T)
D3_barcode$cell_barcode <- str_split_fixed(D3_barcode$cell_barcode,pattern = '_',n = 2)[,2]
D4_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD4.csv',header = T)
D4_barcode$cell_barcode <- str_split_fixed(D4_barcode$cell_barcode,pattern = '_',n = 2)[,2]

D5_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD5.csv',header = T)
D5_barcode$cell_barcode <- str_split_fixed(D5_barcode$cell_barcode,pattern = '_',n = 2)[,2]
D6_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD6.csv',header = T)
D6_barcode$cell_barcode <- str_split_fixed(D6_barcode$cell_barcode,pattern = '_',n = 2)[,2]

D7_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD7.csv',header = T)
D7_barcode$cell_barcode <- str_split_fixed(D7_barcode$cell_barcode,pattern = '_',n = 2)[,2]
D8_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD8.csv',header = T)
D8_barcode$cell_barcode <- str_split_fixed(D8_barcode$cell_barcode,pattern = '_',n = 2)[,2]

D9_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD9.csv',header = T)
D9_barcode$cell_barcode <- str_split_fixed(D9_barcode$cell_barcode,pattern = '_',n = 2)[,2]
D10_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD10.csv',header = T)
D10_barcode$cell_barcode <- str_split_fixed(D10_barcode$cell_barcode,pattern = '_',n = 2)[,2]

D11_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD11.csv',header = T)
D11_barcode$cell_barcode <- str_split_fixed(D11_barcode$cell_barcode,pattern = '_',n = 2)[,2]
D12_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD12.csv',header = T)
D12_barcode$cell_barcode <- str_split_fixed(D12_barcode$cell_barcode,pattern = '_',n = 2)[,2]

D16_barcode <- read.csv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/csv/rawdataD16.csv',header = T)
D16_barcode$cell_barcode <- str_split_fixed(D16_barcode$cell_barcode,pattern = '_',n = 2)[,2]


merged_barcode <- rbind(D1_barcode,D2_barcode,D3_barcode,D4_barcode,D5_barcode,D6_barcode,
                        D7_barcode,D8_barcode, D9_barcode,D10_barcode,D12_barcode,D11_barcode)








####--tsv--
D1_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D1clusters.tsv', col_names = T) %>% as.data.frame() #
D2_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D2clusters.tsv', col_names = T) %>% as.data.frame() #
D2_tsv <- D2_tsv[D2_tsv$assignment %in% c(0,1),]
rownames(D2_tsv) <- 1:nrow(D2_tsv)

D3_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D3clusters.tsv', col_names = T) %>% as.data.frame() #
D4_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D4clusters.tsv', col_names = T) %>% as.data.frame() #
D4_tsv <- D4_tsv[D4_tsv$assignment %in% c(0,1),]
rownames(D4_tsv) <- 1:nrow(D4_tsv)

D5_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D5clusters.tsv', col_names = T) %>% as.data.frame() #
D6_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D6clusters.tsv', col_names = T) %>% as.data.frame() #
D6_tsv <- D6_tsv[D6_tsv$assignment %in% c(0,1),]
rownames(D6_tsv) <- 1:nrow(D6_tsv)

D7_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D7clusters.tsv', col_names = T) %>% as.data.frame() #
D8_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D8clusters.tsv', col_names = T) %>% as.data.frame() #
D8_tsv <- D8_tsv[D8_tsv$assignment %in% c(0,1),]
rownames(D8_tsv) <- 1:nrow(D8_tsv)

D9_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D9clusters.tsv', col_names = T) %>% as.data.frame() #
D10_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D10clusters.tsv', col_names = T) %>% as.data.frame() #
D10_tsv <- D10_tsv[D10_tsv$assignment %in% c(0,1),]
rownames(D10_tsv) <- 1:nrow(D10_tsv)

D11_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D11clusters.tsv', col_names = T) %>% as.data.frame() #
D12_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D12clusters.tsv', col_names = T) %>% as.data.frame() #
D12_tsv <- D12_tsv[D12_tsv$assignment %in% c(0,1),]
rownames(D12_tsv) <- 1:nrow(D12_tsv)

D11_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D11clusters.tsv', col_names = T) %>% as.data.frame() #
D12_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D12clusters.tsv', col_names = T) %>% as.data.frame() #
D12_tsv <- D12_tsv[D12_tsv$assignment %in% c(0,1),]
rownames(D12_tsv) <- 1:nrow(D12_tsv)


D15_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D15clusters.tsv', col_names = T) %>% as.data.frame() #
D16_tsv <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/tsv/D16clusters.tsv', col_names = T) %>% as.data.frame() #
D16_tsv <- D16_tsv[D16_tsv$assignment %in% c(0,1),]
rownames(D16_tsv) <- 1:nrow(D16_tsv)



D1_2clusters <- read_tsv(file = 'D:/IFLT-CLT/2.Analysis/integrated/SNP/D1_D2/clusters.tsv', col_names = T) %>% as.data.frame()



####---sampleD1-2----
EP_barcode <- D1_barcode[D1_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]
PR_barcode <- D2_barcode[D2_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]


EP_TCR <- D1_TCR[D1_TCR$barcode %in% EP_barcode$cell_barcode, ] ## 
PR_TCR <- D2_TCR[D2_TCR$barcode %in% PR_barcode$cell_barcode, ] ## 

PR_tsv <- D2_tsv

PR_0_tsv <- subset(PR_tsv, assignment == 0) 
PR_1_tsv <- subset(PR_tsv, assignment == 1)

PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]

 #Comparing double stranded merged genes
T_cell_clone_uniq <- EP_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
EP_TCR <- customer_clone


T_cell_clone_uniq <- PR_0_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <-  dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_0_TCR <- customer_clone



T_cell_clone_uniq <- PR_1_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <-  dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(5) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_1_TCR <- customer_clone




####---sampleD3-4----
EP_barcode <- D3_barcode[D3_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]
PR_barcode <- D4_barcode[D4_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]

EP_TCR <- D3_TCR
EP_TCR <- D3_TCR[D3_TCR$barcode %in% EP_barcode$cell_barcode, ] 
PR_TCR <- D4_TCR[D4_TCR$barcode %in% PR_barcode$cell_barcode, ] 

PR_tsv <- D4_tsv
PR_0_tsv <- subset(PR_tsv, assignment == 0) 
PR_1_tsv <- subset(PR_tsv, assignment == 1)

PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]

#Comparing double stranded merged genes
T_cell_clone_uniq <- EP_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
EP_TCR <- customer_clone


T_cell_clone_uniq <- PR_0_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_0_TCR <- customer_clone



T_cell_clone_uniq <- PR_1_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(5) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_1_TCR <- customer_clone



####---sampleD5-6----
 EP_barcode <- D5_barcode[D5_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                      "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                      "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                      "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                      "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                      "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                      "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                      "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                      "NK-c5-MCM high", "NK-c6-MT high"),]
 PR_barcode <- D6_barcode[D6_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                      "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                      "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                      "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                      "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                      "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                      "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                      "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                      "NK-c5-MCM high", "NK-c6-MT high"),]
 
 EP_TCR <- D5_TCR
 EP_TCR <- D5_TCR[D5_TCR$barcode %in% EP_barcode$cell_barcode, ] 
 PR_TCR <- D6_TCR[D6_TCR$barcode %in% PR_barcode$cell_barcode, ] 
 
 PR_tsv <- D6_tsv
 
 PR_0_tsv <- subset(PR_tsv, assignment == 0) 
 PR_1_tsv <- subset(PR_tsv, assignment == 1)
 
 PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
 PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]
 
 #Comparing double stranded merged genes
 T_cell_clone_uniq <- EP_TCR
 multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
 T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
 T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
 TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
 library(parallel)
 cl <- makeCluster(3) 
 barcodes_cells <- TCR_info$barcode %>% unique
 clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
 clusterEvalQ(cl, library(dplyr))
 customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
   tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
   clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
   # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
   dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
 }) %>% do.call(rbind, .)
 stopCluster(cl)
 EP_TCR <- customer_clone
 
 
 T_cell_clone_uniq <- PR_0_TCR
 multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
 T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
 T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
 TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
 library(parallel)
 cl <- makeCluster(3) 
 barcodes_cells <- TCR_info$barcode %>% unique
 clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
 clusterEvalQ(cl, library(dplyr))
 customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
   tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
   clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
   # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
   dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
 }) %>% do.call(rbind, .)
 stopCluster(cl)
 PR_0_TCR <- customer_clone
 
 
 
 T_cell_clone_uniq <- PR_1_TCR
 multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
 T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
 T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
 TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
 library(parallel)
 cl <- makeCluster(5) 
 barcodes_cells <- TCR_info$barcode %>% unique
 clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
 clusterEvalQ(cl, library(dplyr))
 customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
   tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
   clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
   # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
   dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
 }) %>% do.call(rbind, .)
 stopCluster(cl)
 PR_1_TCR <- customer_clone
 
 

####---sampleD7-8----
EP_barcode <- D7_barcode[D7_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]
PR_barcode <- D8_barcode[D8_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]

EP_TCR <- D7_TCR
EP_TCR <- D7_TCR[D7_TCR$barcode %in% EP_barcode$cell_barcode, ] 
PR_TCR <- D8_TCR[D8_TCR$barcode %in% PR_barcode$cell_barcode, ] 

PR_tsv <- D8_tsv

PR_0_tsv <- subset(PR_tsv, assignment == 0) 
PR_1_tsv <- subset(PR_tsv, assignment == 1)

PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]

#Comparing double stranded merged genes
T_cell_clone_uniq <- EP_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)
}) %>% do.call(rbind, .)
stopCluster(cl)
EP_TCR <- customer_clone


T_cell_clone_uniq <- PR_0_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_0_TCR <- customer_clone



T_cell_clone_uniq <- PR_1_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(5) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_1_TCR <- customer_clone




####---sampleD9-10----
EP_barcode <- D9_barcode[D9_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]
PR_barcode <- D10_barcode[D10_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]

EP_TCR <- D9_TCR
EP_TCR <- D9_TCR[D9_TCR$barcode %in% EP_barcode$cell_barcode, ] ## 
PR_TCR <- D10_TCR[D10_TCR$barcode %in% PR_barcode$cell_barcode, ] 

PR_tsv <- D10_tsv

PR_0_tsv <- subset(PR_tsv, assignment == 0) 
PR_1_tsv <- subset(PR_tsv, assignment == 1)

PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]


#Comparing double stranded merged genes
T_cell_clone_uniq <- EP_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)###
}) %>% do.call(rbind, .)
stopCluster(cl)
EP_TCR <- customer_clone


T_cell_clone_uniq <- PR_0_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <-  dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_0_TCR <- customer_clone



T_cell_clone_uniq <- PR_1_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] 
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <-  dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(5) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_1_TCR <- customer_clone










####---sampleD11-12----
EP_barcode <- D11_barcode[D11_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                     "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                     "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                     "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                     "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                     "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                     "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                     "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                     "NK-c5-MCM high", "NK-c6-MT high"),]
PR_barcode <- D12_barcode[D12_barcode$cell_type %in% c("CD4-c1-SELL", "CD4-c2-CCR6", 
                                                       "CD4-c3-GZMK", "CD4-c4-XCL1", "CD4-c5-GNLY", "CD4-c6-FOXP3", 
                                                       "CD4-c7-MT2A", "CD4-c8-IFRD1", "CD8-c1-SELL", "CD8-c2-IL7R", 
                                                       "CD8-c3-GZMK", "CD8-c4-FGFBP2", "CD8-c5-TYROBP", "CD8-c6-LNPEP", 
                                                       "CD8-c7-IFIT2", "CD8-c8-HSPA6", 
                                                       "dgT-c1-TRDV1", "dgT-c2-TRDV2", "dgT-c3-LNPEP",
                                                       "MAIT-c1", "MAIT-c2-MT high", "MAIT-c3-HSP high",
                                                       "NK-c1-FGFBP2", "NK-c2-GZMK", "NK-c3-SELL", "NK-c4-LNPEP", 
                                                       "NK-c5-MCM high", "NK-c6-MT high"),]

EP_TCR <- D11_TCR
EP_TCR <- D11_TCR[D11_TCR$barcode %in% EP_barcode$cell_barcode, ] 
PR_TCR <- D12_TCR[D12_TCR$barcode %in% PR_barcode$cell_barcode, ] 

PR_tsv <- D12_tsv

PR_0_tsv <- subset(PR_tsv, assignment == 0) 
PR_1_tsv <- subset(PR_tsv, assignment == 1)

PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]


#Comparing double stranded merged genes
T_cell_clone_uniq <- EP_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
EP_TCR <- customer_clone


T_cell_clone_uniq <- PR_0_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)### 
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_0_TCR <- customer_clone



T_cell_clone_uniq <- PR_1_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(5) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_1_TCR <- customer_clone


####---sampleD14-16----
EP_barcode <- D15_tsv$barcode
PR_barcode <- D16_barcode

EP_TCR <- D15_TCR
EP_TCR <- D15_TCR[D15_TCR$barcode %in% EP_barcode, ] 
PR_TCR <- D16_TCR[D16_TCR$barcode %in% PR_barcode$cell_barcode, ] 
PR_tsv <- D16_tsv

PR_0_tsv <- subset(PR_tsv, assignment == 0) 
PR_1_tsv <- subset(PR_tsv, assignment == 1)

PR_0_TCR <- PR_TCR[PR_TCR$barcode %in% PR_0_tsv$barcode,] 
PR_1_TCR <- PR_TCR[PR_TCR$barcode %in% PR_1_tsv$barcode,]


T_cell_clone_uniq <- EP_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(, 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  # dat <- data.frame(cellbarcode = x, clonetype_customer = clonetype_customer)
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)
}) %>% do.call(rbind, .)
stopCluster(cl)
EP_TCR <- customer_clone


T_cell_clone_uniq <- PR_0_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(., 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(3) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_0_TCR <- customer_clone



T_cell_clone_uniq <- PR_1_TCR
multiple_chains <- T_cell_clone_uniq[duplicated(T_cell_clone_uniq$barcode), ] %>% `[`(., 1) %>% unique()
T_cell_clone_uniq <- T_cell_clone_uniq[T_cell_clone_uniq$barcode %in% multiple_chains, ] ##301031 qualited T cells
T_cell_clone_uniq$unique_clone <-  with(T_cell_clone_uniq, paste0(v_gene, j_gene, cdr3))
TCR_info <- dplyr::select(T_cell_clone_uniq, c(barcode, chain, unique_clone)) 
library(parallel)
cl <- makeCluster(5) 
barcodes_cells <- TCR_info$barcode %>% unique
clusterExport(cl = cl, varlist = c("TCR_info", "barcodes_cells"))
clusterEvalQ(cl, library(dplyr))
customer_clone <- pblapply(cl = cl, X = barcodes_cells, FUN = function(x){
  tmp <- TCR_info %>% subset(barcode == x) %>% select(chain, unique_clone) %>% arrange(desc(chain), desc(unique_clone))
  clonetype_customer <- tmp$unique_clone %>% paste0(., collapse = "")
  dat <- c(cellbarcode = x, clonetype_customer = clonetype_customer)
}) %>% do.call(rbind, .)
stopCluster(cl)
PR_1_TCR <- customer_clone






