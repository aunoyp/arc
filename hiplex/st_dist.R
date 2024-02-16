## ----eval=FALSE-------------------------------------------------------------------------------------------------------------
## current_file <- rstudioapi::getActiveDocumentContext()$path
## output_file <- stringr::str_replace(current_file, '.Rmd', '.R')
## knitr::purl(current_file, output = output_file)
## file.edit(output_file)


## ---------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(tictoc)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)
library(png)
library(cowplot)
library(magick)
library(scales)
library(packcircles)
library(ggalt)
library(cccd)
library(stringr)
library(ggfortify)
library(circlize)
library(reshape2)
library(igraph)
library(ggpubr)
library(ggstatsplot)
library(rstatix)
#library("leiden")
source("/home/aunoy/st/arc_profiling/st_analysis/src/st_functions.R", echo=FALSE)


## ---------------------------------------------------------------------------------------------------------------------------
data_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/rethresholded'
clump_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/clumps'
meta_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/overlay'
output_dir_plot = '/home/aunoy/st/arc_profiling/st_analysis/results/plots'
output_dir_tbls = '/home/aunoy/st/arc_profiling/st_analysis/results/tables'


## ---------------------------------------------------------------------------------------------------------------------------
meta_ntrscts = read.csv(file.path(clump_dir, 'meta', 'META_ntrsct.csv'), header = FALSE) %>%
  as_tibble()


## ---------------------------------------------------------------------------------------------------------------------------
df_408 <- load_slice(slice = "408")
df_408 <- clean_408(df_408)
df_408 <- embed_horizontal_408(df_408)
rownames(df_408) = 1:nrow(df_408)
jy_408 = df_408 %>%
  dplyr::select(-c(area, IMAGE.NAME, X_horz, Y_horz, clump)) %>%
  t() %>%
  CreateSeuratObject()
jy_408$image <- df_408$IMAGE.NAME
jy_408 <- NormalizeData(jy_408, scale.factor = 1e5) ###
normed = GetAssayData(jy_408, slot = 'data')
normed[normed < 3] = 0
#for(gene_name in rownames(jy_408)) {
#  mdn_gene_expr = median(normed[gene_name, normed[gene_name, ] > 0])
#}
jy_408 <- SetAssayData(jy_408, slot = 'data', normed)
jy_408 <- FindVariableFeatures(jy_408, selection.method = "vst")
all.genes <- rownames(jy_408)
jy_408 <- ScaleData(jy_408, features = all.genes)
jy_408 <- RunPCA(jy_408, approx = FALSE)
jy_408 <- FindNeighbors(jy_408, dims = 1:30)
jy_408 <- FindClusters(jy_408, resolution = 0.8)
jy_408 <- RunUMAP(jy_408, dims = 1:30)
jy_408$clump = df_408$clump
hcoords = df_408 %>% dplyr::select(c('X_horz', 'Y_horz')) %>% as.matrix()
colnames(hcoords) <- c('pixel_1', 'pixel_2')

jy_408[["H"]] <- CreateDimReducObject(embeddings = hcoords, key = "pixel_", assay = DefaultAssay(jy_408))

df_164 <- load_slice(slice = "164")
df_164 <- clean_164(df_164)
df_164 <- embed_horizontal_164(df_164)

rownames(df_164) = 1:nrow(df_164)
jy_164 = df_164 %>%
  dplyr::select(-c(area, IMAGE.NAME, X, Y, clump)) %>%
  t() %>%
  CreateSeuratObject()
jy_164$image = df_164$IMAGE.NAME

jy_164 <- NormalizeData(jy_164, scale.factor = 1e5) ###
normed = GetAssayData(jy_164, slot = 'data')
normed[normed < 3] = 0
#for(gene_name in rownames(jy_164)) {
#  mdn_gene_expr = median(normed[gene_name, normed[gene_name, ] > 0])
#}
jy_164 <- SetAssayData(jy_164, slot = 'data', normed)

jy_164 <- FindVariableFeatures(jy_164, selection.method = "vst")
all.genes <- rownames(jy_164)
jy_164 <- ScaleData(jy_164, features = all.genes)
jy_164 <- RunPCA(jy_164, approx = FALSE)
jy_164 <- FindNeighbors(jy_164, dims = 1:30)
jy_164 <- FindClusters(jy_164, resolution = 0.8)
jy_164 <- RunUMAP(jy_164, dims = 1:30)
jy_164$clump <- df_164$clump

unique(df_164$IMAGE.NAME)
images_ordered = c('TC_Cortical3', 'TC_Cortical2', 'TC_Cortical1', 'TC_10', 'TC_9', 'TC_8', 'TC_7', 'TC_6', 'TC_5', 'TC_4', 'TC_3', 'TC_2','TC_1','CC_2','CC_3',
           'CC_4', 'CC_5', 'CC_6', 'CC_7', 'CC_9', 'CC_8', 'CC_10',
           'CC_L2-1', 'CC_L2-2', 'CC_L2-3', 'CC_Cortical1', 'CC_Cortical2')

x_horz = 1:length(images_ordered) * 35
y_horz = rep(0, length(images_ordered))
horz_embedding = data.frame()
df_164$X_horz = -1
df_164$Y_horz = -1

images = list.files(meta_dir)
for(i in 1:length(images_ordered)){
    image_name = images_ordered[i]
    print(image_name)
    split_names = strsplit(image_name, '_')
    cortex = toupper(split_names[[1]][1])
    number = split_names[[1]][2]
    number_csv = paste0('_', number, '.csv')
    filename = images[grepl(cortex, images) & grepl(number_csv, images) & grepl('164', images)]
    coordinates = read.table(file.path(meta_dir, filename), sep = ',', header = TRUE)
    if(image_name == "CC_L2-1"){
      coordinates = coordinates[c(1:37, 39:nrow(coordinates)), ]
    }
    ## checked already that lists are equal, missing 1, 18, 19 for now, layer 1 and others
 
    ## so this is a little tricky, so need to get it right
    ## Remember, it is the top right that the coordinate is coming from, but
    ## the bottom right is the new coordinate space.
    ## so first when we get the original coordinate space, to set to relative
    ## of bottom would be the same X, but 1024 - Y
    
    ## push out the coordinates for better visualization
    #x_repelled <- (512 - coordinates$X_Coordinate_In_pixels)
    
    
    df_164[df_164$IMAGE.NAME == image_name, 'X_horz'] = (coordinates$X_Coordinate_In_pixels / 
                                                      IMAGE_SIZE * IMAGE_LEN) + y_horz[i]
    df_164[df_164$IMAGE.NAME == image_name, 'Y_horz'] = ((1024-coordinates$Y_Coordinate_In_pixels) / 
                                                      IMAGE_SIZE * IMAGE_LEN) + x_horz[i]
}

hcoords = df_164 %>% dplyr::select(c('X_horz', 'Y_horz')) %>% as.matrix()
colnames(hcoords) <- c('pixel_1', 'pixel_2')

jy_164[["H"]] <- CreateDimReducObject(embeddings = hcoords, key = "pixel_", assay = DefaultAssay(jy_164))

jy_164<- RenameCells(jy_164, c(outer('164_', 1:ncol(jy_164), FUN=paste0)))
jy_164$area = df_164$area
jy_408<- RenameCells(jy_408, c(outer('408_', 1:ncol(jy_408), FUN=paste0)))
jy_408$area = df_408$area
jy_all <- merge(jy_164, jy_408)

jy_all <- NormalizeData(jy_all, scale.factor = 1e5) ###
normed = GetAssayData(jy_all, slot = 'data')
normed[normed < 3] = 0
for(gene_name in rownames(jy_all)) {
  if (gene_name == 'DCX'){
    mdn_gene_expr = 0.5
    print('skip dcx')
  } else if (!gene_name %in% c('COUPTF2', 'SP8')){
    mdn_gene_expr = median(normed[gene_name, normed[gene_name, ] > 0])
  }else{
    mdn_gene_expr = quantile(normed[gene_name, normed[gene_name, ] > 0], .40)
  }

  normed[gene_name, normed[gene_name, ] < mdn_gene_expr] = 0
}
jy_all <- SetAssayData(jy_all, slot = 'data', normed)

jy_all <- FindVariableFeatures(jy_all, selection.method = "vst")
all.genes <- rownames(jy_all)
jy_all <- ScaleData(jy_all, features = all.genes)
jy_all <- RunPCA(jy_all, approx = FALSE)
jy_all <- FindNeighbors(jy_all, dims = 1:30)
jy_all <- FindClusters(jy_all, resolution = 0.5)
jy_all <- RunUMAP(jy_all, dims = 1:30)

new.cluster.ids = c('CGE/LGE',
                    'Ex',
                    'TBR1+ CGE',
                    'CALB2+DLX2+',
                    'VIP+GAD1+',
                    'SST+LHX6+',
                    'MGE')

names(new.cluster.ids) <- levels(jy_all)
jy_all <- RenameIdents(jy_all, new.cluster.ids)


## ---------------------------------------------------------------------------------------------------------------------------
## Kind of simple
## For each image, make a graph, and then ha


## ---------------------------------------------------------------------------------------------------------------------------
jy_all$slice <- 'err'

jy_all$slice[1:ncol(jy_164)] <- '164'
jy_all$slice[(ncol(jy_164)+1):ncol(jy_all)] <- '408'

jy_all$image_slice <- paste0(jy_all$slice, '_', jy_all$image)


## ---------------------------------------------------------------------------------------------------------------------------
h_all <- rbind(Embeddings(jy_164, "H"), Embeddings(jy_408, "H"))
rownames(h_all) <- colnames(jy_all)

all_embed <-CreateDimReducObject(embeddings = h_all, key = "pixel_", assay = DefaultAssay(jy_all))
jy_all[["H"]] <- all_embed


## ---------------------------------------------------------------------------------------------------------------------------
jy_all <- ProjectDim(jy_all, reduction = "H")

## ---------------------------------------------------------------------------------------------------------------------------
fimage <- jy_all$image_slice[1]
test_image <- jy_all[, which(jy_all$image_slice == fimage)]


## ---------------------------------------------------------------------------------------------------------------------------
test_embed <- Embeddings(test_image, 'H')
d_mat <- dist(test_embed, diag = TRUE, upper = TRUE)
net <- graph_from_adjacency_matrix(as.matrix(d_mat), mode = "undirected", weighted = TRUE)

## ---------------------------------------------------------------------------------------------------------------------------
net <- simplify(net, remove.multiple = F, remove.loops = T)
plot(net)

## ---------------------------------------------------------------------------------------------------------------------------
hist(E(net)$weight)


## ---------------------------------------------------------------------------------------------------------------------------
dmat_df <- as.data.frame(cbind(E(net)$weight, fimage))
colnames(dmat_df) <- c("dist", "image")
dmat_df$dist <- as.numeric(dmat_df$dist)


## ---------------------------------------------------------------------------------------------------------------------------
dmat_df %>%
  ggplot(aes(x = dist)) + geom_histogram()

## ---------------------------------------------------------------------------------------------------------------------------
dist_df <- get_dist_df(jy_all)


## ----fig.width = 15, fig.height = 3-----------------------------------------------------------------------------------------
dist_df %>%
  ggplot(aes(x = image, y = dist)) + geom_violin() + RotatedAxis()

## ---------------------------------------------------------------------------------------------------------------------------
vms_labels <- c(paste0("164_", antr_TC_MS), paste0("408_", post_TC_MS))
dms_labels <- c(paste0("164_", antr_CC_MS), paste0("408_", post_CC_MS))
cc_labels <- c(paste0("164_", antr_CC_Cx), paste0("408_", post_CC_Cx))
tc_labels <- c(paste0("164_", antr_TC_Cx), paste0("408_", post_TC_Cx))


## ---------------------------------------------------------------------------------------------------------------------------
dist_df$region <- 'error'
dist_df$region[which(dist_df$image %in% vms_labels)] = 'VMS'
dist_df$region[which(dist_df$image %in% dms_labels)] = 'DMS'
dist_df$region[which(dist_df$image %in% tc_labels)] = 'VCx'
dist_df$region[which(dist_df$image %in% cc_labels)] = 'DCx'


## ---------------------------------------------------------------------------------------------------------------------------
my_comparisons <- list( c("VMS", "DMS"))
dist_df %>%
  ggplot(aes(x = region, y = dist, color = region)) + geom_boxplot() + RotatedAxis() +   stat_compare_means(comparisons = my_comparisons)


## ---------------------------------------------------------------------------------------------------------------------------
quantile(dist_df$dist)


## ---------------------------------------------------------------------------------------------------------------------------
### Instead of a df, I think it might be better to get a graph here
g <- get_dist_graph(jy_all, verbose = FALSE)
V(g)$image_slice <- jy_all$image_slice

## ---------------------------------------------------------------------------------------------------------------------------
E(g)[from("164_64_64")]$weight

## ---------------------------------------------------------------------------------------------------------------------------
g1 <- get_dist_graph(jy_all, "164_CC_8")

## ---------------------------------------------------------------------------------------------------------------------------
dd <- E(g1)[from("164_37_37")]


## ---------------------------------------------------------------------------------------------------------------------------
get_dist_interactions(g1, jy_all, classes, c_df)


## ---------------------------------------------------------------------------------------------------------------------------

get_dist_interactions <- function(g, obj, classes, c_df, px_thr = 7){
    c_df$val = 0
    for(i in 1:length(V(g))){
      cell_name = V(g)[i]$name
      class_cell = Idents(obj)[which(colnames(obj) == cell_name)]
      ## Get the cells that are the same
      edge_obj <- E(g)[from(cell_name)]
      edge_dists <- edge_obj$weight
      itrct_ix <- edge_dists <= px_thr
      itrcts = as.data.frame(ends(g, edge_obj, 
                                  names = TRUE))[['V2']][which(itrct_ix)]
      neighbor_idents = Idents(obj)[which(colnames(obj) %in% itrcts)]
      for(neighbor in neighbor_idents){
        value = c_df$val[which(c_df$c1 == class_cell & c_df$c2 == neighbor)]
        c_df$val[which(c_df$c1 == class_cell & c_df$c2 == neighbor)] = value + 1
      }
    }
    return(c_df)
}
### Iterate through all images
### run 100 permutations, and get the p value from the real data
permutations = 99
neighbors = 5
threshold = 4

## IMages
imgs <- unique(jy_all$image_slice)

### Need to set up the dataframe
classes <- levels(Idents(jy_all))
n = length(classes)
c1 = rep(classes, n)
c2 = rep(classes, each=n)
c_df <- as.data.frame(cbind(c1, c2))

## will save all into 2 different matrices
interaction_mat <- matrix(-1, nrow(c_df), length(imgs)) 
avoidance_mat <- matrix(-1, nrow(c_df), length(imgs)) 

## Need to create a jy_all embedding
embed_164 <- Embeddings(jy_164, "H")
embed_408 <- Embeddings(jy_408, "H")
jy_all[["H"]] <- CreateDimReducObject(embeddings = rbind(embed_164, embed_408), key = "pixel_", assay = DefaultAssay(jy_all))

for(i in 1:length(imgs)){
  print(paste("Permuting image", i, "of", length(imgs)))
  img <- imgs[i]
  obj <- jy_all[, which(jy_all$image_slice == img)]
  ### Create nn graph for real
  interaction_graph <- get_dist_graph(obj, img)
  ### Iterate through each cell type
  ### Compute the mean for this one
  data_cts <- get_dist_interactions(interaction_graph, obj, classes, c_df, threshold)
  obj_idents <- Idents(obj)
  perm_matrix = matrix(NA, nrow(c_df), permutations)
  for(permutation in 1:permutations){
    Idents(obj) <- sample(Idents(obj))
    perm_cts <- get_dist_interactions(interaction_graph, obj, classes, c_df)
    perm_matrix[, permutation] <- perm_cts$val
  }
  interaction_mat[, i] <- rowSums(perm_matrix <= data_cts$val) / (permutations + 1)
  avoidance_mat[, i] <- rowSums(perm_matrix >= data_cts$val) / (permutations + 1)
  ## count the identities
  tbl <- table(Idents(obj))
  ## If the cell type does not exist in the image, do not count the P value for the
  ## interaction... obviously
  nonexist_ixs <- which(!(c_df$c1 %in% names(tbl)) | !(c_df$c2 %in% names(tbl)))
  interaction_mat[nonexist_ixs, i] <- NA
  avoidance_mat[nonexist_ixs, i] <- NA
}

## ---------------------------------------------------------------------------------------------------------------------------
## Save this difficultly made files
saveRDS(interaction_mat, file.path('/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/', 'interaction_dmat_4px.rds'))
saveRDS(avoidance_mat, file.path('/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/', 'avoidance_dmat_4px.rds'))


## ---------------------------------------------------------------------------------------------------------------------------
interaction_mat <- readRDS(file.path('/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/', 'interaction_dmat_4px.rds'))
avoidance_mat <- readRDS(file.path('/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/', 'avoidance_dmat_4px.rds'))


## ---------------------------------------------------------------------------------------------------------------------------
rownames(interaction_mat) <- paste0(c_df$c1, "->", c_df$c2)
rownames(avoidance_mat) <- paste0(c_df$c1, "->", c_df$c2)
colnames(interaction_mat) <- imgs
colnames(avoidance_mat) <- imgs


## ---------------------------------------------------------------------------------------------------------------------------
plot_map <- interaction_mat
plot_map[which(plot_map >= 0.05)] = NA
pheatmap(plot_map, cluster_cols = FALSE, cluster_rows = FALSE, cex = 0.8,
         cell_width = 5, cell_height = 5)


## ---------------------------------------------------------------------------------------------------------------------------
sum(interaction_mat < 0.05, na.rm = TRUE)


## ---------------------------------------------------------------------------------------------------------------------------
min(rowMeans(interaction_mat, na.rm = TRUE))


## ---------------------------------------------------------------------------------------------------------------------------
rownames(avoidance_mat) <- paste0(c_df$c1, "->", c_df$c2)
colnames(avoidance_mat) <- imgs
rownames(interaction_mat) <- paste0(c_df$c1, "->", c_df$c2)
colnames(interaction_mat) <- imgs


## ---------------------------------------------------------------------------------------------------------------------------
layouts = c(layout_with_dh,
            layout_with_drl,
            layout_nicely,
            layout_with_fr,
            layout_with_graphopt,
            layout_with_lgl,
            layout_with_gem
            )
layout_names = c('dh', 'drl', 'nicely', 'fr', 'graphopt', 'lgl', 'gem')


## ---- fig.height = 15, fig.width = 15---------------------------------------------------------------------------------------
plot_ig <- function(g, lt_fxn, lt_name, stream){
  return(plot(g, layout=lt_fxn(net), main= paste(stream, lt_name)))
}

vms_labels <- c(paste0("164_", antr_TC_MS), paste0("408_", post_TC_MS))
dms_labels <- c(paste0("164_", antr_CC_MS), paste0("408_", post_CC_MS))

streams = c('vms', 'dms')

for(stream in streams){

  if(stream == 'vms'){
    stream_labels = vms_labels
  } else{
    stream_labels = dms_labels
  }

  sub_mat <- interaction_mat[, stream_labels]
  avgs <- rowMeans(sub_mat < 0.05, na.rm = TRUE)
  
  interaction_df <- c_df %>% add_column(value = avgs)
  colnames(interaction_df) <- c('from', 'to', 'value')
  mat <- dcast(interaction_df, from ~ to, value.var = "value") %>%
    column_to_rownames("from") %>%
    #mutate_all(.funs = readr::parse_number) %>%
    as.matrix() 
  net <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)
  
  #### Get correct positions
  class_ixs <- c()
  for(class in V(net)$name ){
    class_ixs <- c(class_ixs, which(classes == class))
  }
  
  V(net)$color <- get_cluster_colors(V(net)$name)
  
  net <- simplify(net, remove.loops = TRUE)
  
  if(all(colnames(sub_mat) == vms_labels)){
    size_scores = vms_scores
  } else{ size_scores = dms_scores}
  
  sw <- c()
  for(vertex in V(net)$name){
    wt <- E(net)[from(vertex) & to(vertex)]$weight
    if(length(wt) == 0){
      sw <- c(sw, 0)
    } else{
    sw <- c(sw, E(net)[from(vertex) & to(vertex)]$weight)
    }
  }
  ## First one sets up sizes as interactions, but that's less important
  #V(net)$size <- 20 + size_scores[,V(net)$name] / 12
  V(net)$size <- avgs[seq(1, nrow(interaction_mat), 8)][class_ixs]*70 + 12
  ### So I need to get the self weight essentially
  
  # The labels are currently node IDs.
  # Setting them to NA will render no labels:
  
  ### First one sets the label size based on self interaction, now it's based on interaction
  #V(net)$label.cex <- avgs[seq(1, nrow(interaction_mat), 8)][class_ixs] + 1 ## Flipped
  V(net)$label.cex <- size_scores[,V(net)$name] / 600 * 2 + 1
  V(net)$label.dist <- 3
  V(net)$frame.width = 2
  V(net)$frame.color = "black"
  V(net)$label.family = 'sans'
  V(net)$label <- NA
  # Set edge width based on weight:
  ## Filter out edges
  net <- delete.edges(net,E(net)[which(E(net)$weight < 0.3)])
  E(net)$width <- 2
  palf <- colorRamp(c(rgb(1,1,1, .2),rgb(0,0,0, .7)), alpha=TRUE)
  edge_cols <- palf(E(net)$weight)
  edge_cols[which(is.nan(edge_cols)[, 1]), ] = 0
  E(net)$color <- rgb(edge_cols, maxColorValue = 255)
  E(net)$lty <- E(net)$weight
  
  
  #change arrow size and edge color:
  E(net)$arrow.size <- 0.5
  #E(net)$edge.color <- "Red"
  if(stream == "vms"){
    vnet = net
  }
  # We can even set the network layout:
  #coords <- layout_(net, nicely())
  #plot(net, layout=layout.star(net, center = V(net)[2]), main="star")
  #par(mfrow=c(3,3))
  #lapply(1:length(layouts), function(i){
  #  plot_ig(net, layouts[[i]], layout_names[[i]], stream)
  #})
}

## ----fig.height = 5, fig.width = 5------------------------------------------------------------------------------------------
#ig.layout <- layout_nicely(vnet)
plot(vnet, layout=ig.layout, main="vms", rescale = TRUE, edge.lty= E(vnet))
rownames(ig.layout) <- V(vnet)$name
plot(net, layout=ig.layout[V(net)$name,], main="dms", edge.lty= E(net), rescale = TRUE)


## ---------------------------------------------------------------------------------------------------------------------------
img = '408_TC_3'
g1 <- get_dist_graph(jy_all, img)
im <- get_dist_interactions(g1, jy_all, classes, c_df, 4)


## ---------------------------------------------------------------------------------------------------------------------------
obj <- jy_all[, which(jy_all$image_slice == img)]
### Create nn graph for real
XY_obj = Embeddings(obj, 'H')
g1=delete.edges(g1, which(E(g1)$weight > 4)) # here's my condition.

## ---------------------------------------------------------------------------------------------------------------------------
lo <- layout.norm(as.matrix(XY_obj))
plot.igraph(g1, 
    layout = lo,
    resize = FALSE) 


## ---------------------------------------------------------------------------------------------------------------------------
vms_mat <- interaction_mat[, vms_labels]
dms_mat <- interaction_mat[, dms_labels]

res <- prop.test(x = c(sum(vms_mat[1, ] < 0.05), sum(dms_mat[1, ] < 0.05)), 
                 n = c(sum(vms_mat[1, ] >= 0.05), sum(dms_mat[1, ] >= 0.05)))
# Printing the results
res 

## ---------------------------------------------------------------------------------------------------------------------------
tbl <- rbind(c(sum(vms_mat[1, ] < 0.05, na.rm = TRUE), sum(dms_mat[1, ] < 0.05, na.rm = TRUE)), 
              c(sum(vms_mat[1, ] >= 0.05, na.rm = TRUE), sum(dms_mat[1, ] >= 0.05, na.rm = TRUE)))
fisher.test(tbl)

## ---------------------------------------------------------------------------------------------------------------------------
rownames(tbl) <- c('hmg', 'htr')
colnames(tbl) <- c('vms', 'dms')

mosaicplot(tbl,
  main = "Mosaic plot",
  color = TRUE
)


## ---------------------------------------------------------------------------------------------------------------------------
#fisher.test

hmg_interaction <- interaction_mat[seq(1, nrow(interaction_mat), 8), ]
htr_interaction <- interaction_mat[!1:nrow(interaction_mat) %in% seq(1, nrow(interaction_mat), 8), ]

tbl <- rbind(c(sum(hmg_interaction < 0.05, na.rm = TRUE), sum(htr_interaction < 0.05, na.rm = TRUE)), 
              c(sum(hmg_interaction >= 0.05, na.rm = TRUE), sum(htr_interaction >= 0.05, na.rm = TRUE)))
fisher.test(tbl)


## ---------------------------------------------------------------------------------------------------------------------------
rownames(tbl) <- c('enriched', 'not-enriched')
colnames(tbl) <- c('self-interaction', 'hetero-interaction')

mosaicplot(tbl,
  main = "Mosaic plot",
  color = TRUE
)

## ---------------------------------------------------------------------------------------------------------------------------
x <- c()
for (row in rownames(tbl)) {
  for (col in colnames(tbl)) {
    x <- rbind(x, matrix(rep(c(row, col), tbl[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
colnames(df) <- c("Status", "Group")
df


## ---------------------------------------------------------------------------------------------------------------------------
test <- fisher.test(table(df))

# combine plot and statistical test with ggbarstats

ggbarstats(
  df, Status, Group,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
)


## ---------------------------------------------------------------------------------------------------------------------------
E()$weight

## ---------------------------------------------------------------------------------------------------------------------------
get_mins <- function(obj, img){
  g <- get_dist_graph(obj, img)
  sobj <- obj[, which(obj$image_slice == img)]
  same_dist = c()
  othr_dist = c()
  cnames <- colnames(sobj)
  for(cname in cnames){
    ### Get sames
    cident = Idents(sobj)[which(colnames(sobj) == cname)]
    sames = cnames[which(Idents(sobj) == cident)]
    E_c <- E(g)[from(cname)]
    if(length(sames) < 2){
      min_same = NA
    } else{
      min_same = min(E_c[to(sames)]$weight)
    }
    othrs = cnames[which(Idents(sobj) != cident)]
    if(length(othrs) == 0){
      min_othrs = NA
    } else{
      min_othrs = min(E_c[to(othrs)]$weight)
    }
    same_dist <- c(same_dist, min_same)
    othr_dist <- c(othr_dist, min_othrs)
  }
  df <- as.data.frame(cbind(cnames, same_dist, othr_dist, as.character(Idents(sobj))))
  colnames(df) <- c('c_id', 'same', 'other', 'ident')
  return(df)
}

## ---------------------------------------------------------------------------------------------------------------------------
tc_3_df <- get_mins(jy_all, '408_TC_3')

## ---------------------------------------------------------------------------------------------------------------------------
na.omit(tc_3_df) %>%
  pivot_longer(!c(c_id, ident), names_to = "type", values_to = "dist") %>%
  ggplot(aes(x = ident, y = as.numeric(dist), color = as.factor(type))) + geom_boxplot()

## ---------------------------------------------------------------------------------------------------------------------------
### Nearest distance df
### https://stackoverflow.com/questions/29402528/append-data-frames-together-in-a-for-loop

datalist = list()
# or pre-allocate for slightly more efficiency
datalist = vector("list", length = length(imgs))

tic()
for (i in 1:length(imgs)) {
    datalist[[i]] <- get_mins(jy_all, imgs[i]) # add it to your list
}
toc()
### Takes 41.756 seconds

nd_df = do.call(rbind, datalist)

## ---------------------------------------------------------------------------------------------------------------------------
nd_df$area <- jy_all$area
nd_df$image <- jy_all$image_slice
nd_df$stream <- 'err'
nd_df$stream[which(nd_df$image %in% dms_labels)] <- 'DMS'
nd_df$stream[which(nd_df$image %in% vms_labels)] <- 'VMS'
nd_df$stream <- as.factor(nd_df$stream)
nd_df$same <- as.numeric(nd_df$same)
nd_df$other <- as.numeric(nd_df$other)


## ---------------------------------------------------------------------------------------------------------------------------
na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident, image, stream), names_to = "type", values_to = "dist") %>%
  filter(stream != 'err') %>%
  ggplot(aes(x = ident, y = as.numeric(dist), color = as.factor(stream))) + geom_boxplot() + facet_wrap(~type) + RotatedAxis()


## ---------------------------------------------------------------------------------------------------------------------------
na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident), names_to = "type", values_to = "dist") %>%
  filter(grepl('MS', area)) %>%
  ggplot(aes(x = area, y = as.numeric(dist), color = as.factor(type))) + geom_violin() + RotatedAxis()

## ---------------------------------------------------------------------------------------------------------------------------
na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident, image, stream), names_to = "type", values_to = "dist") %>%
  filter(stream != 'err') %>%
  ggplot(aes(x = stream, y = as.numeric(dist), color = as.factor(type))) + geom_boxplot() + RotatedAxis()


## ---------------------------------------------------------------------------------------------------------------------------
stat.test <- na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident, image, stream), names_to = "type", values_to = "dist") %>%
  filter(stream != 'err') %>%
  filter(ident == 'CGE/LGE') %>%
  mutate(type = as.factor(type)) %>%
  dplyr::select(c('stream', 'type', 'dist')) %>%
  group_by(stream) %>%
  t_test(dist ~ type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test


## ---------------------------------------------------------------------------------------------------------------------------
df <- na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident, image, stream), names_to = "type", values_to = "dist") %>%
  filter(stream != 'err') %>%
  filter(ident == 'CGE/LGE') %>%
  mutate(type = as.factor(type)) %>%
  dplyr::select(c('stream', 'type', 'dist'))
bxp <- ggboxplot(
  df, x = "stream", y = "dist", 
  color = "type", palette = c("#00AFBB", "#E7B800")
  )

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "stream", dodge = 0.8)
bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  )

# Add 10% spaces between the p-value labels and the plot border
bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


## ---------------------------------------------------------------------------------------------------------------------------
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj", tip.length = 0,
  remove.bracket = TRUE
  )

# Show adjusted p-values and significance levels
# Hide ns (non-significant)
bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0, hide.ns = TRUE
  )

## ---------------------------------------------------------------------------------------------------------------------------
df$stream <- droplevels(df$stream)
stat.test2 <- df %>%
  t_test(dist ~ stream, p.adjust.method = "bonferroni")
stat.test2

## ---------------------------------------------------------------------------------------------------------------------------
stat.test <- stat.test %>%
  add_xy_position(x = "stream", dodge = 0.8)
bxp.complex <- bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  )
# 2. Add stat.test2
# Add more space between brackets using `step.increase` 
stat.test2 <- stat.test2 %>% add_xy_position(x = "stream")
bxp.complex <- bxp.complex + 
  stat_pvalue_manual(
    stat.test2,  label = "p", tip.length = 0.02,
    step.increase = 0.05
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# 3. Display the plot
bxp.complex 

## ----fig.width = 4, fig.height = 2------------------------------------------------------------------------------------------
bxp.complex 


## ---------------------------------------------------------------------------------------------------------------------------
na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident, image, stream), names_to = "type", values_to = "dist") %>%
  filter(stream == 'DMS') %>%
  ggplot(aes(dist, after_stat(density), colour = type)) + geom_freqpoly(binwidth = 1)

## ---------------------------------------------------------------------------------------------------------------------------
na.omit(nd_df) %>%
  pivot_longer(!c(c_id, area, ident, image, stream), names_to = "type", values_to = "dist") %>%
  filter(stream == 'DMS') %>%
  dplyr::select(dist) %>%
  as.matrix() %>%
shapiro.test()

## ---------------------------------------------------------------------------------------------------------------------------
clump_df$clump_type = "NaN"
clump_names = str_remove(clump_df$...1, "_Clump")
clump_df$clump_type[which(clump_names %in% linear)] = "Linear"
clump_df$clump_type[which(clump_names %in% psuedolinear)] = "Psuedo-Linear"
clump_df$clump_type[which(clump_names %in% big_clump)] = "True Clump"


## ---------------------------------------------------------------------------------------------------------------------------
cdf <- clump_df %>% filter(clump_type != "NaN")
cdf2 <- cdf[, 2:13]


## ---------------------------------------------------------------------------------------------------------------------------
clump_pca <- prcomp(cdf2, center = TRUE,scale. = TRUE)

summary(clump_pca)

## ---------------------------------------------------------------------------------------------------------------------------
autoplot(clump_pca, data = cdf, colour = 'clump_type')

## ---------------------------------------------------------------------------------------------------------------------------
clump_df %>%
  ggplot(aes(x = clump_type, y = eval(as.symbol("Perim.")))) + geom_boxplot()


## ---------------------------------------------------------------------------------------------------------------------------
### Time for Ripley's analysis



## ---------------------------------------------------------------------------------------------------------------------------


### Iterate through all images
### run 100 permutations, and get the p value from the real data

### Need to set up the dataframe
classes <- levels(Idents(jy_all))
n = length(classes)
c1 = rep(classes, n)
c2 = rep(classes, each=n)
c_df <- as.data.frame(cbind(c1, c2))

## Need to create a jy_all embedding
embed_164 <- Embeddings(jy_164, "H")
embed_408 <- Embeddings(jy_408, "H")
jy_all[["H"]] <- CreateDimReducObject(embeddings = rbind(embed_164, embed_408), key = "pixel_", assay = DefaultAssay(jy_all))

vms_scores <- matrix(0, 1, length(classes))
dms_scores <- matrix(0, 1, length(classes))
colnames(vms_scores) <- classes
colnames(dms_scores) <- classes

for(i in 1:length(imgs)){
  #print(paste("Permuting image", i, "of", length(imgs)))
  img <- imgs[i]
  if(!is.na(interaction_mat['CGE/LGE->CGE/LGE', img]) && interaction_mat['CGE/LGE->CGE/LGE', img] < 0.01){
    obj <- jy_all[, which(jy_all$image_slice == img)]
    ### Create nn graph for real
    interaction_graph <- get_dist_graph(obj, img)
    ### Iterate through each cell type
    ### Compute the mean for this one
    data_cts <- get_dist_interactions(interaction_graph, obj, classes, c_df, threshold)
    scores <- data_cts$val[seq(1, nrow(interaction_mat), 8)]
    #score <- data_cts$val[1]
    #tbl <- table(Idents(obj))
    #norm_score <- score / tbl['CGE/LGE'] ### This shoudl never be NA
    if(img %in% vms_labels){
      vms_scores <- vms_scores + scores
    }
    if(img %in% dms_labels){
      dms_scores <- dms_scores + scores
    }
  }
}


## ---------------------------------------------------------------------------------------------------------------------------
hist(vms_scores)

## ---------------------------------------------------------------------------------------------------------------------------
hist(dms_scores)

## ---------------------------------------------------------------------------------------------------------------------------
### So if X and Y represent random variables that draw from the same distribution
t.test(vms_scores, y = dms_scores)      # P = .00001855

## ---------------------------------------------------------------------------------------------------------------------------
jy_all.markers <- FindAllMarkers(jy_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
jy_all.markers %>%
   group_by(cluster) %>%
   slice_max(n = 32, order_by = avg_log2FC)

