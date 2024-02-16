## st_profiling.R
## Aunoy Poddar
## 06/08/2022
## Script adapted from st_profiling notebook to generate the plots of interest

## ------------------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(tictoc)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)

#### AUNOY
# I am going to want to add which areas are in teh analysis, how many cells, and date as a meta text file


## ------------------------------------------------------------------------------------------------------------------------------------------
## Define functions to use
### log(x+1) function
log1 <- function(x) {
  return(log(x+1))
}

### Scale and then log+1
log1_and_mult <- function(x) {
  return(log1(x*1e5))
}


## ------------------------------------------------------------------------------------------------------------------------------------------
data_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/rethresholded'
output_plot_dir = '/home/aunoy/st/arc_profiling/st_analysis/results/plots'
output_tble_dir = '/home/aunoy/st/arc_profiling/st_analysis/results/tables'

## ------------------------------------------------------------------------------------------------------------------------------------------

### right now i want experiment # for the day and overwrite to be false

## Write to file
date <- gsub('-', '', Sys.Date())
exp_num <- '2'
folder <- paste0(date, '_', exp_num)

## check if folder exists, if not, then proceed
expltdir <- file.path(output_plot_dir, folder)
extbldir <- file.path(output_tble_dir, folder)

if(dir.exists(expltdir) || dir.exists(extbldir)){
  quit('Folder exists. Set --overwrite to TRUE to replace folder')
}

## Create the directories
dir.create(expltdir)
dir.create(extbldir)

## ------------------------------------------------------------------------------------------------------------------------------------------
## Load data from data_dir
print('Loading data...')
df = data.frame()
for (file_name in list.files(data_dir)){
  df_to_append <- read.table(file.path(data_dir, file_name), sep = ',', header = TRUE)
  df_to_append <- df_to_append %>%
    dplyr::select(-1)
  if(grepl('164', file_name)){
    df_to_append <- df_to_append %>%
      dplyr::select(-X)
  }
  colnames(df_to_append) <- toupper(colnames(df_to_append))
  df_to_append <- df_to_append %>%
    mutate(area = strsplit(file_name, '.csv')[[1]])
  if(!is_empty(df)){
    df_to_append <- df_to_append %>%
      dplyr::select(colnames(df))
  }
  df <- rbind(df, df_to_append)
}

## ------------------------------------------------------------------------------------------------------------------------------------------
print('Writing metainformation to file...')
fileConn<-file(file.path(output_tble_dir, 'meta.txt'), open = "w")
exp_id       <-  paste0('experiment id\t:\t', paste0(date, '-', exp_num))
current_time  <- paste0('script started\t:\t', Sys.time())
cell_num      <- paste0('# of cells used\t:\t', as.character(nrow(df)))
areas_used    <- paste0('areas analyzed\t:\t', paste0(unique(df$area), collapse = ', '))
writeLines(c(exp_id, current_time,cell_num, areas_used), fileConn)
close(fileConn)

## ------------------------------------------------------------------------------------------------------------------------------------------
## Name cells arbitrarily numbers
rownames(df) <- c(outer(c('C'), 1:dim(df)[1], FUN=paste0))

## ------------------------------------------------------------------------------------------------------------------------------------------
## Define gene lists for cell origin / phenomena
migra = toupper(c('Dcx', 'Lrp8', 'Reln', 'Dcdc2', 'Ncam1', 'Kia0319', 'Vldlr'))
#CGE = toupper(c('Egfr', 'Vip', 'Prox1'))
#LGE = toupper(c('Tshz1', 'Gsx2', 'Emx1'))
#MGE = toupper(c('Lhx6', 'Maf1', 'Sst'))
CGE = toupper(c('Egfr', 'Prox1'))
LGE = toupper(c('Tshz1', 'Gsx2', 'Emx1'))
MGE = toupper(c('Lhx6', 'Maf1'))
CGE_LGE = toupper(c('Scgn', 'Couptf2', 'Sp8', 'Calb2', 'Pax6'))
GABA = toupper(c('Dlx2', 'Gad1'))
mature_IN = toupper(c('Gad1', 'Vip', 'Sst'))
progen_IN = toupper(c('Dlx2'))
Excit = toupper(c('Eomes', 'Tbr1', 'Satb2'))
ligand = toupper(c('Reln', 'Cxcl12', 'Cxcl14'))
recept = toupper(c('Lrp8', 'Cxcr7', 'Cxcr4', 'Vldlr'))


## ------------------------------------------------------------------------------------------------------------------------------------------
## Pivot wider and longer for various plotting and analysis
df_longer <- df %>%
  rownames_to_column('Cell') %>%
  dplyr::select(-area) %>%
  pivot_longer(!Cell, names_to = 'Gene', values_to = "Puncta2Nuc_IR")

df_wide <- df %>%
  rownames_to_column('Cell')%>%
  dplyr::select(-area) %>%
  pivot_longer(!Cell, names_to = 'Gene', values_to = "Puncta2Nuc_IR") %>%
  pivot_wider(names_from = Gene, values_from = Puncta2Nuc_IR)

## ------------------------------------------------------------------------------------------------------------------------------------------

print("Plotting heatmaps...")
### Get the numeric values only
df_num <- df_wide %>%
  column_to_rownames('Cell')

### Normalize to cell total
### Log + Multiply
### Center and scale
scaled_mat <- df_num %>%
  sweep(1, rowSums(df_num), '/') %>%
  log1_and_mult() %>%
  apply(1, scale)

### Apply separate clustering, need to look into why 
hr <- hclust(as.dist(1-cor(t(scaled_mat), method="pearson")), method = "complete")
hc <- hclust(as.dist(1-cor(scaled_mat, method="spearman")), method="complete")

annotation <- data.frame(area = factor(df$area))
rownames(annotation) <- rownames(df) # check out the row names of annotation

rownames(scaled_mat) <- colnames(df_num)
scaled_mat %>%
  pheatmap(annotation_col = annotation, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  #breaks = seq(from= -5, to = 5, by = 11/100.), 
  cluster_rows = hr, cluster_cols = hc, 
  #cluster_rows = TRUE, cluster_cols = TRUE, 
  #clustering_distance_rows = "correlation", 
  #clustering_distance_cols = "correlation",
  fontsize_row = 8, fontsize_col = 4, show_colnames = FALSE,
           show_rownames = TRUE, cell_width = 6, cellheight = 8,
           filename = file.path(output_plot_dir, 'pheatmap_by_dendro.png'))


## ------------------------------------------------------------------------------------------------------------------------------------------
rownames(scaled_mat) <- colnames(df_num)
scaled_mat %>%
  pheatmap(annotation_col = annotation, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  #breaks = seq(from= -5, to = 5, by = 11/100.), 
  cluster_rows = hr, cluster_cols = FALSE, 
  #cluster_rows = TRUE, cluster_cols = TRUE, 
  #clustering_distance_rows = "correlation", 
  #clustering_distance_cols = "correlation",
  fontsize_row = 8, fontsize_col = 4, show_colnames = FALSE,
           show_rownames = TRUE, cell_width = 6, cellheight = 8,
           filename = file.path(output_plot_dir, 'pheatmap_by_area.png'))


## ------------------------------------------------------------------------------------------------------------------------------------------
dist_mat <- as.dist(1-cor(scaled_mat, method="spearman"))
dist_mat %>%
pheatmap(annotation_col = annotation,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_rows = TRUE, cluster_cols = TRUE,
fontsize_row = 8, fontsize_col = 4, show_colnames = FALSE,
         show_rownames = FALSE, cell_width = 0.1, cellheight = 0.3,
          filename = file.path(output_plot_dir, 'pheatmap_by_cell2cell_distance.png'))


## ------------------------------------------------------------------------------------------------------------------------------------------
## plot out the histograms
print("Plotting histograms")

plots = list()
all_genes <- df %>% 
  dplyr::select(-area) %>%
  colnames()

for (i in 1:length(all_genes)){
  plots[[i]] <- df_longer %>%
  filter(Gene == all_genes[i]) %>%
  mutate(log_norm = log1_and_mult(Puncta2Nuc_IR)) %>%
  ggplot(aes(x=log_norm, label = Gene)) + 
  geom_histogram(color="black", fill="white") + 
  labs(title=all_genes[i])+
  theme_classic()
}

#gridExtra::grid.arrange(grobs = plots, ncol = 4, nrow = 8, lengths=2:6)
#if (save){
ml <- marrangeGrob(plots, nrow=2, ncol=2)
ggsave(file.path(output_plot_dir, 'lognorm_histograms.pdf'), ml)
#} else{
#  marrangeGrob(plots, nrow=2, ncol=2)
#}

# ------------------------------------------------------------------------------------------------------------------------------------------
print("Generating and processing Seurat object...")
jyobj <- df %>%
  dplyr::select(-area) %>%
  t() %>%
  CreateSeuratObject()


## ------------------------------------------------------------------------------------------------------------------------------------------
jyobj <- NormalizeData(jyobj) ###


## ------------------------------------------------------------------------------------------------------------------------------------------
print("Plotting varianle features...")

jyobj <- FindVariableFeatures(jyobj, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(jyobj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(jyobj) + theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 6))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 6))
plot1 + plot2 + theme(axis.title.y = element_text(size = 10))
ggsave(file.path(output_plot_dir, 'vst_plot.jpg'))

# ## ------------------------------------------------------------------------------------------------------------------------------------------
all.genes <- rownames(jyobj)
jyobj <- ScaleData(jyobj, features = all.genes)
# 
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
jyobj <- RunPCA(jyobj, features = VariableFeatures(object = jyobj), approx = FALSE)
# print(jyobj[["pca"]], dims = 1:5, nfeatures = 5)
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
# VizDimLoadings(jyobj, dims = 1:2, reduction = "pca")
# 
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
jyobj$area <- df$area
# DimPlot(jyobj, reduction = "pca", group.by = 'area')
# 
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
# DimHeatmap(jyobj, dims = 1,  balanced = TRUE)
# 
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
# DimHeatmap(jyobj, dims = 1:15, balanced = TRUE)
# 
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
print("Creating JackStraw Plot...")

jyobj <- JackStraw(jyobj, num.replicate = 100)
jyobj <- ScoreJackStraw(jyobj, dims = 1:20)
# 
JackStrawPlot(jyobj, dims = 1:15)
ggsave(file.path(output_plot_dir, 'jackstraw.png'))
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
# ElbowPlot(jyobj)
# 
# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
# 
## do a grid search on dims to use for clustering
# resolution https://arxiv.org/pdf/0803.0476.pdf

print("Creating UMAP plots through grid search...")

pcdims = c(5, 10, 15, 20, 25, 30)
resolutions = c(0.1, 0.3, 0.5, 0.8, 1, 3, 10)
i = 1

tic()
for (resolution in resolutions){
  for (pcdim in pcdims){
    jyobj <- FindNeighbors(jyobj, dims = 1:pcdim)
    jyobj <- FindClusters(jyobj, resolution = resolution)
    jyobj <- RunUMAP(jyobj, dims = 1:pcdim)
    plots[[i]] = DimPlot(jyobj, reduction = "umap", 
                         group.by = 'seurat_clusters') +
      ggtitle(paste0("PCdim:", as.character(pcdim), "|Res:", as.character(resolution)))
    plots[[i+1]] = DimPlot(jyobj, reduction = "umap", group.by = 'area') +
      ggtitle(paste0("PCdim:", as.character(pcdim), "|Res:", as.character(resolution)))
    i = i + 2
  }
}
ml <- marrangeGrob(plots, nrow=1, ncol=1)
ggsave(file.path(output_plot_dir, 'umap_plots.pdf'), ml)
toc()

## ------------------------------------------------------------------------------------------------------------------------------------------
print("Plotting features...")

plots = list()
for (i in 1:length(all_genes)){
  plots[[i]] <- FeaturePlot(jyobj, pt.size = 0.5,features = all_genes[i])
}

ml <- marrangeGrob(plots, nrow=2, ncol=2)
ggsave(file.path(output_plot_dir, 'feature_plots.pdf'), ml)

# 
# ## ------------------------------------------------------------------------------------------------------------------------------------------
print("Writing differential gene expression table...")

jyobj.markers <- FindAllMarkers(jyobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
jyobj.markers %>%
   group_by(cluster) %>%
   slice_max(n = 32, order_by = avg_log2FC)
write.table(jyobj.markers, file = file.path(output_tble_dir, "diff_gene_expr_table.txt"), sep = "\t")

# ## ------------------------------------------------------------------------------------------------------------------------------------------

print("Plotting coexpression...")

pairs <- data.frame(rbind(c("EGFR", "COUPTF2"),
                          c("PROX1", "SCGN"),
                          c("SP8", "COUPTF2"),
                          c("DCX", "GAD1"),
                          c("TBR1", "SATB2"),
                          c("TBR1", "DCX"),
                          c("RELN", "VLDLR"),
                          c("SCGN", "VLDLR"),
                          c("VLDLR", "COUPTF2"),
                          c("VLDLR", "DCX"),
                          c("TBR1", "GAD1")))

plots = list()
for (i in 1:nrow(pairs)){
  plotobj <- FeaturePlot(jyobj, features = as.character(pairs[i, ]), blend = TRUE)
  idx <- i*4-3
  plots[[idx]] <- plotobj[[1]]
  plots[[idx+1]] <- plotobj[[2]]
  plots[[idx+2]] <- plotobj[[3]]
  plots[[idx+3]] <- plotobj[[4]]
}

ml <- marrangeGrob(plots, nrow=2, ncol=2)
ggsave(file.path(output_plot_dir, "coexpression.pdf"), ml)

# ## ------------------------------------------------------------------------------------------------------------------------------------------
print("Plotting raw vs log...")

plots = list()
for (i in 1:length(all_genes)){
  
  scaled_mat <- df_num %>%
    sweep(1, rowSums(df_num), '/') %>%
    log1_and_mult() %>%
    apply(1, scale)
  rownames(scaled_mat) <- colnames(df_num)
  
  nscale_mat <- df_num %>%
    sweep(1, rowSums(df_num), '/') %>%
    log1_and_mult()
  
  log_t_scale <- scaled_mat[all_genes[i], ]
  log_t_nscale <- nscale_mat[, all_genes[i]]
  raw <- df[, all_genes[i]]
  
  ex_df <- as.data.frame(cbind(log_t_scale,log_t_nscale, raw))
  density_plot = TRUE
  if(i == 19 || i == 23){
    density_plot = FALSE
  }
  
  if(density_plot){
  plots[[i*2-1]] <- ex_df %>%
    ggplot(aes(x = raw, y = log_t_nscale)) + geom_point(size=0.5) + xlim(0, 0.25) + 
    geom_density_2d() + geom_vline(xintercept = 0.005, linetype = 'dashed', 
                                   color = 'red', size = .5)
    ggtitle(paste0(all_genes[i], ' Log No Scale'))
  
  plots[[i*2]] <- ex_df %>%
    ggplot(aes(x = raw, y = log_t_scale)) + geom_point(size=0.5) + xlim(0, 0.25) + 
    geom_density_2d() + geom_vline(xintercept = 0.005, linetype = 'dashed', 
                                   color = 'red', size = .5) +
    ggtitle(paste0(all_genes[i], 'Log Scaled'))
  }
  else{
    plots[[i*2-1]] <- ex_df %>%
      ggplot(aes(x = raw, y = log_t_nscale)) + geom_point(size=0.5) + xlim(0, 0.25) + 
      geom_vline(xintercept = 0.005, linetype = 'dashed', 
                                     color = 'red', size = .5) +
      ggtitle(paste0(all_genes[i], ' Log No Scale'))
    
    plots[[i*2]] <- ex_df %>%
      ggplot(aes(x = raw, y = log_t_scale)) + geom_point(size=0.5) + xlim(0, 0.25) + 
      geom_vline(xintercept = 0.005, linetype = 'dashed', 
                                     color = 'red', size = .5) +
      ggtitle(paste0(all_genes[i], 'Log Scaled'))
  }
}

ml <- marrangeGrob(plots, nrow=1, ncol=2)
ggsave(file.path(output_plot_dir, 'rawvlog.pdf'), ml, width = 6, height = 4.5)

fileConn<-file(file.path(output_tble_dir, 'meta.txt'), open = "a")
endtime = paste0("script finished\t:\t", Sys.time())
writeLines(c(endtime), fileConn)
close(fileConn)