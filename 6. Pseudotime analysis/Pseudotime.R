library(Seurat)
library(Matrix)
library(ggplot2)
library(parallel)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)
library(reshape2)
library(rlang)
library(future)
library(ggsci)
library(SoupX)
library(harmony)
library(patchwork)
library(scales)
library(tidyverse)
library(ggpubr)
library(viridis)
library(clusterProfiler)
library(ghibli)
library(monocle3)
library(reticulate)
library(sceasy)
library("Nebulosa")
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library('GSVA')
library("limma")
library(GSEABase)
library(miloR)
library(SCENIC)
library(AUCell)
library(GENIE3)
library(RcisTarget)


# ----monocle3----

# ----step1 constructing objects => new_cell_data_set()---------
CD8_T_d <- subset(TNK_v7, final_clusters %in% c("CD8-c1-CCR7", "CD8-c2-GZMB", "CD8-c3-GZMK"))

expression_matrix <- GetAssayData(CD8_T_d, assay = 'RNA', slot = 'counts')
cell_metadata <- CD8_T_d@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))

rownames(gene_annotation) <- rownames(expression_matrix)
CD8_T_mono <- new_cell_data_set(expression_matrix,
                              cell_metadata = cell_metadata,
                              gene_metadata = gene_annotation)


#----step2 pre => preprocess_cds-----------
CD8_T_mono <- preprocess_cds(CD8_T_mono, num_dim = 3) 
CD8_T_mono <- align_cds(CD8_T_mono, alignment_group = "orig.ident") 

#----step3 visualize => reduce_dimension()-----------
CD8_T_mono <- reduce_dimension(CD8_T_mono) 

p1 <- plot_cells(CD8_T_mono, label_groups_by_cluster=FALSE,  color_cells_by = "final_clusters")
ggsave('CD8_monocle3_final_clusters_umap.png', plot = p1, width = 6, height = 6, dpi = 200)


#----step4 recluster-----------
CD8_T_mono <- cluster_cells(CD8_T_mono, cluster_method = 'louvain')
p3 <- plot_cells(CD8_T_mono, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p4 <- plot_cells(CD8_T_mono, color_cells_by = "time_point", show_trajectory_graph = FALSE) + ggtitle("label by time_point")
plot_all <- p3+p4
ggsave('CD8_T monocle3_cluster.png', plot = plot_all, width = 12, height = 6, dpi = 200)



#----step6 identify the trajectory => learn_graph()-----------
CD8_T_mono <- learn_graph(CD8_T_mono)


pfind <- plot_cells(CD8_T_mono, cell_size = 0.6,
                    color_cells_by = "final_clusters",
                    label_cell_groups = FALSE, 
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE, 
                    label_branch_points = FALSE) + scale_color_manual(breaks = c("CD8-c1-CCR7", "CD8-c2-GZMB", 
                                                                                 "CD8-c3-GZMK", "CD8-c4-TYMS", "CD8-c5-LAG3"),
                                                                      values = color_used[4:25])


p6 <- plot_cells(CD8_T_mono, cell_size = 0.6,
                 color_cells_by = "final_clusters",
                 label_cell_groups = FALSE, 
                 label_groups_by_cluster = FALSE,
                 label_leaves = FALSE, 
                 label_branch_points = FALSE) + scale_color_manual(breaks = c("CD8-c1-CCR7", "CD8-c2-GZMB", 
                                                                              "CD8-c3-GZMK", "CD8-c4-TYMS", "CD8-c5-LAG3"),
                                                                   values = color_used[4:25])+
  theme(panel.border = element_rect(colour = 'black', size = 1, fill = NA),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  theme(legend.text=element_text(size = 20),
        legend.title=element_text(size = 20, face = 'bold')) + labs(title = '')


ggsave('monocle3_trajectory_find.png', plot = pfind, width = 10, height = 6, dpi = 300)
ggsave('monocle3_trajectory_cluster.png', plot = p6, width = 10, height = 6, dpi = 300)



#----step7 order => order_cells()-----------

p7 <- pfind + geom_hline(yintercept = c(seq(5,5.5,0.1))) + 
  geom_vline(xintercept = c(seq(-9,-8,0.1))) + 
  ggtitle("find root cell")
ggsave('monocle3_root_cell.png', plot = p7, width = 10, height = 6, dpi = 300)


root_cell_screen <- data.frame(CD8_T_mono@int_colData$reducedDims$UMAP)
root_cell_screen1 <- subset(root_cell_screen, X1 > -8.8 & X1 < -8.6 & X2 > 5.2 & X2 < 5.3)

root.cell <- rownames(root_cell_screen1)
CD8_T_mono <- order_cells(CD8_T_mono, root_cells = root.cell)

p8 <- plot_cells(CD8_T_mono,
                 cell_size = 0.6,
                 color_cells_by = "pseudotime", 
                 label_groups_by_cluster = FALSE,
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE,  
                 label_branch_points = FALSE)+ 
  scale_color_viridis(option = "D", begin = 0, end = 1)+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.border = element_rect(colour = 'black', size = 1, fill = NA))+
  labs(color = 'Pseudotime')+
  theme(legend.text=element_text(size = 16),
        legend.title=element_text(size = 18, face = 'bold')) + labs(title = '')

ggsave('monocle3_trajectory_graph_final.png', plot = p8, width = 7, height = 5, dpi = 300)
ggsave('monocle3_trajectory_graph_merge.png', plot = p6+p8, width = 16, height = 5, dpi = 300)



p9 <- plot_cells(CD8_T_mono,
                 cell_size = 0.6,
                 color_cells_by = "Donor", 
                 label_groups_by_cluster = FALSE,
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE,  
                 label_branch_points = FALSE)+ 
  scale_color_viridis(option = "D", begin = 0, end = 1)+
  scale_color_manual(breaks = c("Fetal1", "Fetal2", "Fetal3", "Adult1",
                                "Adult2", "Adult3"),
                     values=c("#07646DFF", "#45837FFF", "#699896FF", "#E394BBFF", 
                              "#DE6DA7FF", "#D4419EFF"))+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.border = element_rect(colour = 'black', size = 1, fill = NA))+
  labs(color = 'Pseudotime')+
  theme(legend.text=element_text(size = 16),
        legend.title=element_text(size = 18, face = 'bold')) + labs(title = '')

ggsave('monocle3_trajectory_graph_donor.png', plot = p9, width = 7, height = 5, dpi = 300)



p10 <- plot_genes_in_pseudotime(CD8_T_mono['ARG1',], color_cells_by="Detail_Annotation_tmp2", 
                                min_expr=0.5, ncol = 1, cell_size = 0.1 )+
  scale_color_manual(breaks = c("CMP/GMP_C28", "Pro neutrophil_C46", "Pre neutrophil_C45",
                                "Neutrophil CAMP_C25", "Neutrophil S100A12_C34", 
                                "Neutrophil ZDHHC19_C10", "Neutrophil FCGR3B_C5"),values=color_tt)
ggsave("Genes_Jitterplot.png", plot = p10, width = 8, height = 4)







