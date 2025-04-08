library(Seurat)
library(SeuratWrappers)
library(tictoc)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
#library(SeuratObject)
library(pracma)
library(monocle3)
library(future)
options(Seurat.Object.assay.version = "v5")
options(future.globals.maxSize = 8000 * 1024^2)

set.seed(12345)

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the "--cores" flag is present
if ("--cores" %in% args) {
  # Get the index of the "--cores" flag
  cores_index <- which(args == "--cores")
  
  # Get the value of the number of cores
  cores <- as.integer(args[cores_index + 1])
} else {
    cores = 1
}

print(cores)
cores = 12
print(cores)

#cds = readRDS("/wynton/group/paredes/Aunoy/integrations_w_IPC/fig2_monocle_cds.rds")
plan("multicore", workers = cores)
plan()

new_partitions <-  function(cds, cell_type_2_remove, threshold = 0.9){
  clusters_to_remove = c()
  bads = unique(cds@clusters$UMAP$clusters[colData(cds)$scvi_class == cell_type_2_remove])
  for(bad in bads) {
      bad_table = table(colData(cds)$scvi_class[cds@clusters$UMAP$clusters == bad])
     percentage = bad_table[cell_type_2_remove] / sum(bad_table)
     print(percentage)
     if(percentage > threshold){
       clusters_to_remove = c(clusters_to_remove, bad)
       print(clusters_to_remove)
     }
  }
  return(clusters_to_remove)
}

# Load interneurons
interneurons = readRDS("/wynton/group/paredes/Aunoy/integrations_w_IPC/arc_shi_tsankova_tsankova_ctx.rds")
interneurons= JoinLayers(interneurons)
interneurons= FindNeighbors(interneurons, reduction = "integrated.scvi", dims = 1:30)
interneurons<- FindClusters(interneurons, resolution = 1.5, cluster.name = "scvi_clusters")
interneurons= RunUMAP(interneurons,reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

inter.exp_cds = SeuratWrappers::as.cell_data_set(interneurons)

# Transfer Seurat UMAP embeddings
inter.exp_cds@int_colData@listData$reducedDims$UMAP <- interneurons@reductions$umap.scvi@cell.embeddings

inter.exp_cds = cluster_cells(reduction_method = "UMAP", cluster_method = "louvain", inter.exp_cds, k = 25)

rec.part <- c(rep(1, length(inter.exp_cds@colData@rownames)))

groups_to_remove = c("LGE D1 MSN PDYN")
partition_num = 2
#for(badgrp in groups_to_remove){
#    bad_clusters = new_partitions(inter.exp_cds, badgrp, threshold = 0.1)
#    rec.part[inter.exp_cds@clusters$UMAP$clusters %in% bad_clusters] = partition_num
#    partition_num = partition_num + 1
#}
#rec.part[inter.exp_cds@clusters$UMAP$clusters %in% c("45", "46")] = 2
names(rec.part) <- inter.exp_cds@colData@rownames
rec.part <- as.factor(rec.part)
inter.exp_cds@clusters$UMAP$partitions <- rec.part
#inter.exp_cds@clusters$UMAP$cluster = colData(inter.exp_cds)$scvi_fate

cds = inter.exp_cds
#remove(inter.exp_cds)


### UNCOMMENT FOR HYPERPARAMETERIZATION
#sus_clusters = cds@clusters$UMAP$cluster

#b <- list()


#nabors = c(30, 40, 50)
#euclids = c(1)
#desics = c(0.4, 0.5, 0.6)
#mbls = c(10, 25, 35)

#count = 1

#for(nabor in nabors) {
#  tic()
#    print(paste("Loop num: ", count))
#    for(mbl in mbls){
            #if(use_clusters){
            #    cds@clusters$UMAP$cluster = colData(cds)$scvi_fate
            #} else{
            #    cds@clusters$UMAP$cluster = sus_clusters
            #}
            
#            cds <- learn_graph(cds, use_partition = FALSE,
#                   close_loop = TRUE, 
#                   learn_graph_control = list(nn.cores = cores, nn.k = nabor,
#                   minimal_branch_len = mbl,
#                euclidean_distance_ratio = 1,
#                geodesic_distance_ratio = 0.4), verbose = TRUE)
#           b[[count]] <-plot_cells(cds,
#           color_cells_by = "cluster",
#           label_groups_by_cluster=FALSE,
#           label_leaves=F,
#           label_branch_points=F,
#           label_roots = F,
#           trajectory_graph_color = "cyan") + 
#           coord_fixed() + ggtitle(paste0("m:", mbl, "|", nabor))
#            count = count + 1
#     #   }
#  }
#  print("timing for graph construction")
#  toc()
#  myplot = cowplot::plot_grid(plotlist = b)
#  ggsave("/wynton/group/paredes/Aunoy/gridsearch_learngraph.pdf", myplot, dpi = 50)
#}

#myplot = cowplot::plot_grid(plotlist = b)
#ggsave("/wynton/group/paredes/Aunoy/gridsearch_learngraph.pdf", myplot)
cds <- learn_graph(cds, use_partition = TRUE, verbose = TRUE, close_loop = TRUE,
                    learn_graph_control = list(nn.cores = cores,
                    nn.k = 40,
                    minimal_branch_len = 35,
                    geodesic_distance_ratio = 0.4))

saveRDS(cds, "/wynton/group/paredes/Aunoy/integrations_w_IPC/fig2_monocle_cds_cleaned.rds")
#
#my.cds = cds
#remove(cds)
#
#my.cds_subset = my.cds[,!is.na(colData(my.cds)$CGE_Cortical)]
#markers <- graph_test(my.cds_subset, neighbor_graph="principal_graph", cores=cores)
#saveRDS(markers, '/wynton/group/paredes/Aunoy/cge_test_graph_test2.rds')
#
#my.cds_subset = my.cds[,!is.na(colData(my.cds)$MGE_Cortical)]
#markers <- graph_test(my.cds_subset, neighbor_graph="principal_graph", cores=cores)
#saveRDS(markers, '/wynton/group/paredes/Aunoy/mge_test_graph_test2.rds')

#my.cds_subset = my.cds[,!is.na(colData(my.cds)$MGE_Striatal_Cells)]
#markers <- graph_test(my.cds_subset, neighbor_graph="principal_graph", cores=cores)
#saveRDS(markers, '/wynton/group/paredes/Aunoy/mge_str_test_graph_test2.rds')

#my.cds_subset = my.cds[,!is.na(colData(my.cds)$LGE_D1_MSN_Lin)]
#markers <- graph_test(my.cds_subset, neighbor_graph="principal_graph", cores=cores)
#saveRDS(markers, '/wynton/group/paredes/Aunoy/lge_d1msn_test_graph_test2.rds')


#my.cds_subset = my.cds[,!is.na(colData(my.cds)$MGE_Subpallial_Cholinergic)]
#markers <- graph_test(my.cds_subset, neighbor_graph="principal_graph", cores=cores)
#saveRDS(markers, '/wynton/group/paredes/Aunoy/mge_SUPCHOL_test_graph_test2.rds')

