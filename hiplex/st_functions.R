### ST Functions
### Aunoy Poddar
### Written Monday Dec 12th to summarize separate functions written throughout
### the course of analysis, which were started in May 2022


### Define the groups that each image is in
post_TC_MS = c('TC_15', 'TC_14', 'TC_13', 'TC_12', 'TC_11', 'TC_10', 'TC_9', 
               'TC_8', 'TC_7', 'TC_6', 'TC_5', 'TC_4', 'TC_3', 'TC_2', 'TC_1')
post_TC_Cx = c('Layer1', 'TC_20', 'TC_19', 'TC_18', 'TC_17', 'TC_16')
post_CC_MS = c('CC_4', 'CC_5', 'CC_6', 'CC_7', 'CC_8', 'CC_9', 'CC_10', 
               'CC_11', 'CC_12')
post_CC_Cx = c('CC_L2-3', 'CC_L2-2', 'CC_L2-1', 'CC_Cortical1', 'CC_Cortical2')
antr_TC_MS = c("TC_1", "TC_2", "TC_3", "TC_4", "TC_5", "TC_6", "TC_7", "TC_8", 
               "TC_9", "TC_10")
antr_TC_Cx = c("TC_Cortical1", "TC_Cortical2", "TC_Cortical3")
antr_CC_MS = c("CC_2", "CC_3", "CC_4", "CC_5","CC_6", "CC_7", "CC_9")
antr_CC_Cx = c("CC_8", "CC_10", "CC_Cortical1", "CC_Cortical2", "CC_L2-1",
               "CC_L2-2", "CC_L2-3")

dms = c(post_CC_MS, antr_CC_MS)
vms = c(post_TC_MS, antr_TC_MS)
cc = c(post_CC_Cx, antr_CC_Cx)
tc = c(post_TC_Cx, antr_TC_Cx)

images_408_ordered = c('TC_20', 'TC_19', 'TC_18', 'TC_17', 'TC_16', 'TC_15', 'TC_14', 'TC_13', 
                   'TC_12', 'TC_11', 'TC_10', 'TC_9', 'TC_8', 'TC_7', 'TC_6', 
                   'TC_5', 'TC_4', 'TC_3', 'TC_2', 'TC_1', 'CC_4', 'CC_5', 
                   'CC_6', 'CC_7', 'CC_8', 'CC_9', 'CC_10', 'CC_11', 
                   'CC_12', 'CC_L2-3', 'CC_L2-2', 'CC_L2-1', 
                   'CC_Cortical1', 'CC_Cortical2')

data_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/rethresholded'
clump_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/clumps'
meta_dir = '/home/aunoy/st/arc_profiling/st_analysis/hand_annotated_data/overlay'

IMAGE_SIZE = 1024
IMAGE_LEN = 25

load_slice <- function(slice, debug = FALSE){
  
  slices = c('408', '164')
  skip = slices[which(slices != slice)]
  df = data.frame()
  
  ### Iterate through the files with the quantifications
  for (file_name in list.files(data_dir)){
    if(grepl(skip, file_name)){
      ### Skip if it is from the anterior slice. Because there are so many 
      next
    }
    if(debug){
      print(file_name)
    }
    df_to_append <- read.table(file.path(data_dir, file_name), 
                               sep = ',', header = TRUE)
    
    while(length(ind <- which(df_to_append$Image.Name %in% c("", "Layer1"))) > 0){
      df_to_append$Image.Name[ind] <- df_to_append$Image.Name[ind -1]
    }
    
    colnames(df_to_append) <- toupper(colnames(df_to_append))
    df_to_append <- df_to_append %>%
      mutate(area = strsplit(file_name, '.csv')[[1]])
    
    ## Add relative_XY_position
    
    if(!is_empty(df)){
      df_to_append <- df_to_append %>%
        dplyr::select(colnames(df))
    }
    df <- rbind(df, df_to_append)
  }
  return(df)
}

clean_408 <- function(df_408, debug = FALSE){
  df_408$IMAGE.NAME = unlist(lapply(df_408$IMAGE.NAME, gsub, 
                                    pattern='_Cluster', replacement=''))
  df_408$IMAGE.NAME = unlist(lapply(df_408$IMAGE.NAME, gsub, 
                                    pattern='[*]', replacement=''))
  df_408$IMAGE.NAME = unlist(lapply(df_408$IMAGE.NAME, gsub, 
                                    pattern='X', replacement=''))
  df_408$IMAGE.NAME = unlist(lapply(df_408$IMAGE.NAME, gsub, 
                                    pattern='L2_', replacement='L2-'))
  df_408$IMAGE.NAME = unlist(lapply(df_408$IMAGE.NAME, gsub, 
                                    pattern='-L2', replacement='_L2'))
  df_408$IMAGE.NAME = unlist(lapply(df_408$IMAGE.NAME, gsub, 
                                    pattern='Tc_12', replacement='TC_12'))
  if (debug){
    print(unique(df_408$IMAGE.NAME))
  }
  return(df_408)
}

clean_164 <- function(df_164, debug = FALSE){
  df_164$IMAGE.NAME = unlist(lapply(df_164$IMAGE.NAME, gsub, pattern='CC-', replacement='CC_'))
  df_164$IMAGE.NAME = unlist(lapply(df_164$IMAGE.NAME, gsub, pattern='[*]', replacement=''))
  df_164$IMAGE.NAME = unlist(lapply(df_164$IMAGE.NAME, gsub, pattern='X', replacement=''))
  df_164$IMAGE.NAME = unlist(lapply(df_164$IMAGE.NAME, gsub, pattern='L2', replacement='CC_L2'))
  df_164$IMAGE.NAME = unlist(lapply(df_164$IMAGE.NAME, gsub, pattern='L2_', replacement='L2-'))
  df_164$IMAGE.NAME = unlist(lapply(df_164$IMAGE.NAME, gsub, pattern='-L2', replacement='_L2'))
  tc_cortical_names_bad = df_164[grepl('TC', df_164$area) & grepl('Cortical', df_164$IMAGE.NAME), 'IMAGE.NAME']
  df_164[grepl('TC', df_164$area) & grepl('Cortical', df_164$IMAGE.NAME), 'IMAGE.NAME'] = unlist(lapply(tc_cortical_names_bad, gsub, pattern='Cort', replacement='TC_Cort'))
  cc_cortical_names_bad = df_164[grepl('CC', df_164$area) & grepl('Cortical', df_164$IMAGE.NAME), 'IMAGE.NAME']
  df_164[grepl('CC', df_164$area) & grepl('Cortical', df_164$IMAGE.NAME), 'IMAGE.NAME'] = unlist(lapply(cc_cortical_names_bad, gsub, pattern='Cort', replacement='CC_Cort'))
  if(debug){
   print(unique(df_164$IMAGE.NAME))
  }
  return(df_164)
}

embed_horizontal_408 <- function(df_408, debug = FALSE){
  ### Set the ordering of the images
  images_ordered = images_408_ordered
  
  ### Separate each group by 35 pixels for plotting
  x_horz = 1:length(images_ordered) * 35
  y_horz = rep(0, length(images_ordered))
  horz_embedding = data.frame()
  
  ### Set column to -1 to check for later issues
  df_408$X_horz = -1
  df_408$Y_horz = -1
  
  ### Assign clumps separately later, so set to NaN for now
  df_408$clump = NaN
  clump_header = '408_'
  clump_files = list.files(clump_dir)
  ## This is the size of an image in the global coordinate space
  IMAGE_LEN = 25
  
  images = list.files(meta_dir)
  for(i in 1:length(images_ordered)){
    image_name = images_ordered[i]
    if(debug){
      print(image_name) 
    }
    split_names = strsplit(image_name, '_')
    cortex = toupper(split_names[[1]][1])
    number = split_names[[1]][2]
    number_csv = paste0('_', number, '.csv')
    filename = images[grepl(cortex, images) & grepl(number_csv, images) & grepl('408', images)]
    coordinates = read.table(file.path(meta_dir, filename), sep = ',', header = TRUE)
    ## in TC_1 and TC_18, JYK removed the 7th cell because the background was too high
    if(image_name %in% c("TC_1", "TC_18")){
      coordinates = coordinates[c(1:6, 8:nrow(coordinates)), ]
    }
    
    df_408[df_408$IMAGE.NAME == image_name, 'X_horz'] = (coordinates$X_Coordinate_In_pixels / 
                                                           IMAGE_SIZE * IMAGE_LEN) + y_horz[i]
    df_408[df_408$IMAGE.NAME == image_name, 'Y_horz'] = ((1024-coordinates$Y_Coordinate_In_pixels) / 
                                                           IMAGE_SIZE * IMAGE_LEN) + x_horz[i]
  }
  return(df_408)
}


embed_horizontal_164 <- function(df_164, debug = FALSE){
  image_names = unique(df_164$IMAGE.NAME)
  # Preset these variables to negative values so I can easily check if they were updated later
  df_164$X = -1
  df_164$Y = -1
  df_164$clump = NaN
  clump_header = '164_'
  clump_files = list.files(clump_dir)
  # set some normalization variables
  ## This is the size of the image when the pixel values are taken from top left down
  IMAGE_SIZE = 1024
  ## This is the size of an image in the global coordinate space
  IMAGE_LEN = 32
  TC_IMAGE_HEIGHT = 410
  TC_IMAGE_WIDTH = 446
  CC_IMAGE_HEIGHT= 422
  CC_IMAGE_WIDTH = 214
  
  # Load the dataframe with global and relative coordinates
  img_cords = read.table(file.path(meta_dir, '164_pixel_coordinates.csv'), sep = ',', header = TRUE)
  
  images = list.files(meta_dir)
  for(image_name in image_names){
    if(grepl('408', image_name)){
      next
    }
    if(debug){
      print(image_name)
    }
    split_names = strsplit(image_name, '_')
    cortex = toupper(split_names[[1]][1])
    number = split_names[[1]][2]
    number_csv = paste0('_', number, '.csv')
    filename = images[grepl(cortex, images) & grepl(number_csv, images) & grepl('164', images)]
    coordinates = read.table(file.path(meta_dir, filename), sep = ',', header = TRUE)
    
    ## Get the clumps
    filename = clump_files[grepl(clump_header, clump_files) & grepl(paste0(image_name, '_'),clump_files)]
    if(length(filename != 0)){
      clump_df = as.data.frame(t(read.csv(file.path(clump_dir, filename), header = FALSE)))
      colnames(clump_df) = c('roi', 'cluster')
      clump_df$roi = as.numeric(clump_df$roi)
      
      if(image_name == "CC_L2-1"){
        coordinates = coordinates[c(1:37, 39:nrow(coordinates)), ]
        if(38 %in% clump_df$roi){clump_df = clump_df[clump_df$roi != 38, ]}
        clump_df$roi[clump_df$roi > 37] = clump_df$roi[clump_df$roi > 37] - 1
      }
      
      image_idxs = which(df_164$IMAGE.NAME == image_name)
      clump_df$roi_idxs = image_idxs[1] + clump_df$roi - 1
      df_164[clump_df$roi_idxs, "clump"] = paste0(clump_header, image_name, '_', clump_df$cluster)
    } else{
      if(debug){
        print('No clumps!')
        print(image_name)
      }
    }
    
    if(cortex == 'CC'){ 
      if(debug){print(paste('cc', filename, image_name))}
      ## So if CC, we add the coordinates for TC_1 to overall image coordinates
      x_adj = img_cords[img_cords$Name == 'TC_1', 'x'] + 
        img_cords[img_cords$Name == 'G_CC1_to_TC1', 'x']
      ## Start from bottom, add the height, subtract TC_1 height, and then global CC1 to TC1
      y_adj = TC_IMAGE_HEIGHT - img_cords[img_cords$Name == 'TC_1', 'y'] +
        img_cords[img_cords$Name == 'G_CC1_to_TC1', 'y'] + CC_IMAGE_HEIGHT
    }else{
      if(debug){print(paste('tc', filename, image_name))}
      x_adj = 0
      y_adj = TC_IMAGE_HEIGHT
    }
    
    ## So don't do repelled for now
    #x_repelled <- (512 - coordinates$X_Coordinate_In_pixels)
    
    ## so the resized x distance is from left, so just add to the box location and adj
    df_164[df_164$IMAGE.NAME == image_name, 'X'] = (coordinates$X_Coordinate_In_pixels / 
                                                      IMAGE_SIZE * IMAGE_LEN) + 
      img_cords[img_cords$Name == image_name, 'x'] + x_adj
    ## resized y distance
    df_164[df_164$IMAGE.NAME == image_name, 'Y'] = y_adj - img_cords[img_cords$Name == image_name, 'y'] - 
      (coordinates$Y_Coordinate_In_pixels / IMAGE_SIZE * 
         IMAGE_LEN)  
  }
  return(df_164)
}

get_slice <- function(sobj){
  return(ifelse(grepl('164', sobj$area[1]), '164', '408'))
}

### Plotting Functions

### This 
plot_clusters_umap <- function(sobj, clustering, pt.size = 3, space = "umap")
{
  coordinates <- Embeddings(sobj, reduction = space)
  expmat  = as.character(Idents(sobj))
  gene_df <- as.data.frame(cbind(coordinates, expmat))
  colnames(gene_df) <- c('X', 'Y', 'expr')
  gene_df$X = as.numeric(gene_df$X)
  gene_df$Y = as.numeric(gene_df$Y)
  summary_gene_df = gene_df %>% dplyr::group_by(expr) %>% dplyr::summarise(xmean = mean(X), ymean = mean(Y))
  plot <- ggplot(gene_df, aes(x = X, y = Y, color = as.factor(expr))) + geom_point(size = pt.size, alpha = 0.8) + 
    theme_classic() + ggtitle(clustering) + NoAxes() + #NoLegend()  + 
    theme(title = element_text(face = 'bold', size = rel(1), hjust = 1)) 
  cluster_colors = scales::hue_pal()(length(unique(expmat)))
  plot = plot + scale_colour_manual(values = cluster_colors)
  return(plot)
}


plot_features_umap <- function(sobj, gene, pt.size = 3, alpha = 0.8, space = "umap", color = '#CB2A55', flipped = FALSE)
{
  coordinates <- Embeddings(sobj, reduction = space)
  expmat <- as.matrix(FetchData(sobj, gene))
  gene_df <- as.data.frame(cbind(coordinates, expmat))
  colnames(gene_df) <- c('X', 'Y', 'expr')
  gene_df <- gene_df %>% dplyr::arrange(!is.na(expr), expr)
  colors = c('grey90', 'grey90', color)
  gene_df$expr[gene_df$expr == 0] = NA
  plot <- ggplot(gene_df, aes(x = X, y = Y, color = expr)) + geom_point(size = pt.size, alpha = alpha)+  
    theme_classic() + ggtitle(gene) +  scale_color_gradient(na.value = colors[1], low = colors[2], high = colors[3], labels = NULL)  + theme(title = element_text(face = 'bold', size = rel(1), hjust = 1)) 
  if(flipped){
    plot <- plot + scale_x_reverse() + scale_y_reverse()
  }
  return(plot + NoAxes() + NoLegend())
}

plot_clusters_vertical_spatial <- function(sobj, cluster, clustering = NULL, anterior = FALSE, 
                                           cluster_color =  '#CB2A55', pt.size = 1, space = "H", 
                                           arc = TRUE, cortical = TRUE)
{
  cluster_identity = Idents(sobj) == cluster
  idents = levels(Idents(sobj))
  coordinates <- Embeddings(sobj, reduction = space)
  gene_df <- as.data.frame(cbind(coordinates, cluster_identity))
  colnames(gene_df) <- c('X', 'Y', 'clust')
  if(arc){
    dms_ixs = which(sobj$image %in% dms | sobj$image %in% cc)
    gene_df$Y[dms_ixs] <- gene_df$Y[dms_ixs] + 50
  }
  
  if(cortical){
    if(get_slice(sobj) == '164'){
      tc = antr_TC_Cx
      cc = antr_CC_Cx
    } else{
      tc = post_TC_Cx
      cc = post_CC_Cx
    }
    non_tc_ix = which(!sobj$image %in% tc)
    cc_ix = which(sobj$image %in% cc)
    gene_df$Y[non_tc_ix] <- gene_df$Y[non_tc_ix] + 15
    gene_df$Y[cc_ix] <- gene_df$Y[cc_ix] + 15
    #min_vms <- min(gene_df$Y[which(sobj$image %in% vms)])
    #max_tc <- max(gene_df$Y[which(sobj$image %in% tc)])
    #max_dms <- max(gene_df$Y[which(sobj$image %in% dms)])
    #min_cc <- min(gene_df$Y[which(sobj$image %in% cc)])
    #boundary_dorsal <- mean(min_vms, max_tc)
    #boundary_ventral <- mean(max_dms, min_cc)
  }
  
  gene_df <- gene_df %>% dplyr::arrange(clust)
  plot <- ggplot(gene_df, aes(x = X, y = Y, color = factor(clust))) + geom_point(size = pt.size, alpha = 1) +  
    theme_classic() + #ggtitle(cluster) + 
    NoAxes() + NoLegend() + 
    coord_fixed(ratio = 0.5)  + theme(title = element_text(face = 'bold', size = rel(0.8), hjust = 1)) 
  cluster_color = scales::hue_pal()(length(idents))[which(idents == cluster)]
  plot = plot + scale_colour_manual(values = c('grey90', cluster_color))
  #if(cortical){plot = plot + 
  #  geom_hline(yintercept=boundary_dorsal, linetype = "dashed",color = cluster_color) +
  #  geom_hline(yintercept=boundary_ventral, linetype = "dashed",color = cluster_color)}
  return(plot)
}

plot_clusters_vertical_spatial_no_grid <- function(sobj, clustering = NULL, anterior = FALSE, 
                                           cluster_color =  '#CB2A55', pt.size = 1, space = "H", 
                                           arc = TRUE, cortical = TRUE, x_width = 35, force_idents = NA)
{
  final_df <- data.frame()
  x_adj = 5
  idents = ifelse(is.na(force_idents),levels(Idents(sobj)), force_idents)
  #x_width = ifelse(get_slice(sobj) == '164', )
  for(cluster in idents){
    cluster_identity = Idents(sobj) == cluster
    coordinates <- Embeddings(sobj, reduction = space)
    gene_df <- as.data.frame(cbind(coordinates, cluster_identity))
    colnames(gene_df) <- c('X', 'Y', 'clust')
    if(arc){
      dms_ixs = which(sobj$image %in% dms | sobj$image %in% cc)
      gene_df$Y[dms_ixs] <- gene_df$Y[dms_ixs] + 50
    }
    
    if(cortical){
      if(get_slice(sobj) == '164'){
        tc = antr_TC_Cx
        cc = antr_CC_Cx
      } else{
        tc = post_TC_Cx
        cc = post_CC_Cx
      }
      non_tc_ix = which(!sobj$image %in% tc)
      cc_ix = which(sobj$image %in% cc)
      gene_df$Y[non_tc_ix] <- gene_df$Y[non_tc_ix] + 15
      gene_df$Y[cc_ix] <- gene_df$Y[cc_ix] + 15
      gene_df$cluster = cluster
      #min_vms <- min(gene_df$Y[which(sobj$image %in% vms)])
      #max_tc <- max(gene_df$Y[which(sobj$image %in% tc)])
      #max_dms <- max(gene_df$Y[which(sobj$image %in% dms)])
      #min_cc <- min(gene_df$Y[which(sobj$image %in% cc)])
      #boundary_dorsal <- mean(min_vms, max_tc)
      #boundary_ventral <- mean(max_dms, min_cc)
    }
    
    gene_df$X = gene_df$X + x_adj
    x_adj = x_adj + x_width
    final_df <- rbind(final_df, gene_df)
  }
  gene_df <- final_df
  gene_df <- gene_df %>% dplyr::arrange(clust)
  gene_df$cluster[which(gene_df$clust == 0)] = NA
  plot <- ggplot(gene_df, aes(x = X, y = Y, color = factor(cluster, levels = idents))) + geom_point(size = pt.size, alpha = 1) +  
    theme_classic() + #ggtitle(cluster) + 
    NoAxes() + NoLegend() + 
    coord_fixed(ratio = 0.5)  + theme(title = element_text(face = 'bold', size = rel(0.8), hjust = 1)) 
  #cluster_color = scales::hue_pal()(length(idents))[which(idents == cluster)]
  if(length(levels(sobj)) != length(idents)){
    idents = idents[which(idents %in% levels(sobj))]
  }
  plot = plot + scale_colour_manual(values = get_cluster_colors(idents), na.value = 'grey90')
  #if(cortical){plot = plot + 
  #  geom_hline(yintercept=boundary_dorsal, linetype = "dashed",color = cluster_color) +
  #  geom_hline(yintercept=boundary_ventral, linetype = "dashed",color = cluster_color)}
  return(plot)
}

plot_clump_celltype <- function(sobj, clump_name,  pt_size = 8) {
  
  sbset_obj = sobj[, which(sobj$image == clump_name)]
  xy = Embeddings(sbset_obj, reduction = 'H')
  expmat  = as.character(Idents(sbset_obj))
  #expmat = FetchData(sbset_obj, gene)
  df <- as.data.frame(cbind(xy, expmat))
  colnames(df) <- c('x', 'y', 'ident')
  colors = hue_pal()(length(levels(Idents(jy_all))))
  colors = colors[which(levels(Idents(jy_all)) %in% df$ident)]
  #colors = c('grey90', 'grey90', color_df$color[color_df$gene == gene])
  df$x = as.numeric(df$x)
  df$y = as.numeric(df$y)
  df$clump = sbset_obj$clump
  # myplot <- 
  df %>%
    ggplot(aes(x = x, y = y, color = factor(ident, levels = levels(Idents(sobj))))) +
    geom_point(size = pt_size) + theme_classic() +# ggtitle(clump_name) + 
    coord_fixed(ratio = 1) + NoAxes() + scale_color_manual(values = colors, name = 'cluster') + 
    geom_encircle(data = filter(df, clump != "NaN"), expand = 0.03,aes(group = clump) ) + xlim(0, 32) + ylim(525, 558) +   theme_cowplot() + NoAxes() + NoLegend()#+ 
  #theme(legend.key.size=unit(0.1,'mm'),
  #    legend.text=element_text(size=4),
  #    legend.title=element_text(size=6))#+ scale_color_manual(values = colors, name = 'cluster', guide = guide_legend(override.aes = list(size=4))) +
  #+ scale_color_gradient(na.value = colors[1], low = colors[2], high = colors[3])
  #+ scale_color_gradient(limits = c(0,11))
  #return(myplot)
}


get_cluster_colors <- function(idents, ncols = 7){
  colors <- hue_pal()(ncols)
  return(colors[match(idents, levels(Idents(jy_all)))])
}


plot_clump_celltype <- function(sobj, clump_name,  pt_size = 8) {
  sbset_obj = sobj[, which(sobj$image == clump_name)]
  xy = Embeddings(sbset_obj, reduction = 'H')
  expmat  = as.character(Idents(sbset_obj))
  #expmat = FetchData(sbset_obj, gene)
  df <- as.data.frame(cbind(xy, expmat))
  colnames(df) <- c('x', 'y', 'ident')
  colors = hue_pal()(length(levels(Idents(jy_all))))
  colors = colors[which(levels(Idents(jy_all)) %in% df$ident)]
  #colors = c('grey90', 'grey90', color_df$color[color_df$gene == gene])
  df$x = as.numeric(df$x)
  df$y = as.numeric(df$y)
  df$clump = sbset_obj$clump
  # myplot <- 
  df %>%
    ggplot(aes(x = x, y = y, color = factor(ident, levels = levels(Idents(sobj))))) +
    geom_point(size = pt_size) + theme_classic() +# ggtitle(clump_name) + 
    coord_fixed(ratio = 1) + NoAxes() + scale_color_manual(values = colors, name = 'cluster') + 
    geom_encircle(data = filter(df, clump != "NaN"), expand = 0.03,aes(group = clump) ) + xlim(0, 32) + ylim(525, 558) + theme_cowplot() + NoAxes() + NoLegend()#+ 
  #theme(legend.key.size=unit(0.1,'mm'),
  #    legend.text=element_text(size=4),
  #    legend.title=element_text(size=6))#+ scale_color_manual(values = colors, name = 'cluster', guide = guide_legend(override.aes = list(size=4))) +
  #+ scale_color_gradient(na.value = colors[1], low = colors[2], high = colors[3])
  #+ scale_color_gradient(limits = c(0,11))
  #return(myplot)
}

save_cluster_image_as_csv <- function(sobj, img_name){
  sub_obj <- sobj[, which(sobj$image == img_name)]
  df <- as.data.frame(cbind(1:ncol(sub_obj), get_cluster_colors(Idents(sub_obj))))
  colnames(df) <- c('roi_name', 'color')
  write_csv(df, file.path(output_dir_tbls, '20221213_1', paste0(img_name, '.csv')))
}

plot_features_vertical_spatial_no_grid <- function(sobj, features, clustering = NULL, anterior = FALSE, 
                                                   color =  '#CB2A55', pt.size = 1, space = "H", alpha = 0.8,
                                                   arc = TRUE, cortical = TRUE, x_width = 35, ft.size = 1.5)
{
  final_df <- data.frame()
  label_df <- data.frame()
  x_adj = 5
  #x_width = ifelse(get_slice(sobj) == '164', )
  for(feature in features){
    expmat <- as.matrix(FetchData(sobj, feature))
    coordinates <- Embeddings(sobj, reduction = space)
    gene_df <- as.data.frame(cbind(coordinates, expmat))
    colnames(gene_df) <- c('X', 'Y', 'expr')
    if(arc){
      dms_ixs = which(sobj$image %in% dms | sobj$image %in% cc)
      gene_df$Y[dms_ixs] <- gene_df$Y[dms_ixs] + 50
    }
    
    if(cortical){
      if(get_slice(sobj) == '164'){
        tc = antr_TC_Cx
        cc = antr_CC_Cx
      } else{
        tc = post_TC_Cx
        cc = post_CC_Cx
      }
      non_tc_ix = which(!sobj$image %in% tc)
      cc_ix = which(sobj$image %in% cc)
      gene_df$Y[non_tc_ix] <- gene_df$Y[non_tc_ix] + 15
      gene_df$Y[cc_ix] <- gene_df$Y[cc_ix] + 15
    }
    
    gene_df$X = gene_df$X + x_adj
    x_adj = x_adj + x_width
    final_df <- rbind(final_df, gene_df)
    lb_df <- gene_df[which.max(gene_df$Y), ]
    lb_df$expr <- feature
    label_df <- rbind(label_df, lb_df)
  }
  gene_df <- final_df
  gene_df <- gene_df %>% dplyr::arrange(!is.na(expr), expr)
  colors = c('grey90', 'grey90', color)
  gene_df$gene[gene_df$expr == 0] = NA
  plot <- ggplot(gene_df, aes(x = X, y = Y, color = expr)) + geom_point(size = pt.size, alpha = alpha) +
    geom_text(data = label_df, x=label_df$X, y=label_df$Y + 35, label=label_df$expr, color = 'black', size = ft.size) + 
    theme_classic() + #ggtitle(cluster) + 
    NoAxes() + NoLegend() + 
    coord_fixed(ratio = 0.5)  + theme(title = element_text(face = 'bold', size = rel(0.8), hjust = 1))
  plot = plot + scale_color_gradient(na.value = colors[1], low = colors[2], high = colors[3])#, labels = NULL)
  return(plot)
}

get_dists <- function(jy_obj){
  ### This function will return from a set of cells, the distances between
  ### all of the cells
  print('nullboy')
}

get_dist_df <- function(jy_obj, px_thr = 6){
  image_slices = unique(jy_obj$image_slice)
  df <- data.frame()
  for(image in image_slices){
    ### get cells in image
    cells <- jy_obj[, which(jy_obj$image_slice == image)]
    embed <- Embeddings(cells, 'H')
    d_mat <- dist(embed, diag = TRUE, upper = TRUE)
    net <- graph_from_adjacency_matrix(as.matrix(d_mat), mode = "undirected", weighted = TRUE)
    dmat_df <- as.data.frame(cbind(E(net)$weight, image))
    colnames(dmat_df) <- c("dist", "image")
    dmat_df$dist <- as.numeric(dmat_df$dist)
    df <- rbind(df, dmat_df)
  }
  return(df)
}


get_dist_graphs <- function(jy_obj, verbose = FALSE){
  image_slices = unique(jy_obj$image_slice)
  for(i in 1:length(image_slices)){
    image = image_slices[i]
    if(verbose){
      print(paste0("Loading image ", i, " of ", length(image_slices),
                   " slices..."))
    }
    ### get cells in image
    cells <- jy_obj[, which(jy_obj$image_slice == image)]
    embed <- Embeddings(cells, 'H')
    d_mat <- dist(embed, diag = TRUE, upper = TRUE)
    net <- graph_from_adjacency_matrix(as.matrix(d_mat), mode = "undirected", weighted = TRUE)
    if(i == 1){
      g <- net
    } else{
      g <- union(g, net, byname = TRUE)
    }
    i = i + 1
  }
  return(g)
}

get_dist_graph <- function(jy_obj, image, verbose = FALSE){
    ### get cells in image
    cells <- jy_obj[, which(jy_obj$image_slice == image)]
    embed <- Embeddings(cells, 'H')
    d_mat <- dist(embed, diag = TRUE, upper = TRUE)
    net <- graph_from_adjacency_matrix(as.matrix(d_mat), mode = "directed", weighted = TRUE)
    return(net)
}

plot_clusters_vertical_spatial_merge <- function(sobj, cluster, stream, clustering = NULL, anterior = FALSE, 
                                           cluster_color =  '#CB2A55', pt.size = 1, space = "H", 
                                           arc = TRUE, cortical = TRUE)
{
  
  if(stream == 'CC'){
    ixs = grepl("CC", sobj$area)
  } else if (stream == 'PirC'){
    ixs = grepl("TC", sobj$area) & grepl("164", sobj$area)
  } else if (stream == 'TC'){
    ixs = grepl("TC", sobj$area) & grepl("408", sobj$area)
  } else{
    stop("Invalid stream name provided")
  }
  
  sobj = sobj[, ixs]
  cluster_identity = Idents(sobj) == cluster
  idents = levels(Idents(sobj))
  coordinates <- Embeddings(sobj, reduction = space)
  gene_df <- as.data.frame(cbind(coordinates, cluster_identity))
  colnames(gene_df) <- c('X', 'Y', 'clust')
  
  if(arc){
    dms_ixs = which(sobj$image %in% dms | sobj$image %in% cc)
    gene_df$Y[dms_ixs] <- gene_df$Y[dms_ixs] + 50
  }
  
  if(cortical){
    if(get_slice(sobj) == '164'){
      tc = antr_TC_Cx
      cc = antr_CC_Cx
    } else{
      tc = post_TC_Cx
      cc = post_CC_Cx
    }
    non_tc_ix = which(!sobj$image %in% tc)
    cc_ix = which(sobj$image %in% cc)
    gene_df$Y[non_tc_ix] <- gene_df$Y[non_tc_ix] + 15
    gene_df$Y[cc_ix] <- gene_df$Y[cc_ix] + 15
    #min_vms <- min(gene_df$Y[which(sobj$image %in% vms)])
    #max_tc <- max(gene_df$Y[which(sobj$image %in% tc)])
    #max_dms <- max(gene_df$Y[which(sobj$image %in% dms)])
    #min_cc <- min(gene_df$Y[which(sobj$image %in% cc)])
    #boundary_dorsal <- mean(min_vms, max_tc)
    #boundary_ventral <- mean(max_dms, min_cc)
  }
  
  ix_164 = startsWith(rownames(gene_df), '164')
  ix_408 = startsWith(rownames(gene_df), '408')
  
  min_164 = min(gene_df$Y[ix_164])
  min_408 = min(gene_df$Y[ix_408])
  
  max_164 = max(gene_df$Y[ix_164])
  max_408 = max(gene_df$Y[ix_408])
  max_avg = (max_164 + max_408) / 5
  
  gene_df$Y[ix_164] = gene_df$Y[ix_164] - min_164
  gene_df$Y[ix_408] = gene_df$Y[ix_408] - min_408
  
  gene_df$Y[ix_164] = (gene_df$Y[ix_164] / (max_164 - min_164)) * max_avg
  gene_df$Y[ix_408] = (gene_df$Y[ix_408] / (max_408 - min_408)) * max_avg
  #gene_df$Y[ix_408] = 0
  
  gene_df <- gene_df %>% dplyr::arrange(clust)
  plot <- ggplot(gene_df, aes(x = X, y = Y, color = factor(clust))) + geom_point(size = pt.size, alpha = 1) +  
    theme_classic() + #ggtitle(cluster) + 
    NoAxes() + NoLegend() + 
    coord_fixed(ratio = 0.5)  + theme(title = element_text(face = 'bold', size = rel(0.8), hjust = 1)) 
  cluster_color = scales::hue_pal()(length(idents))[which(idents == cluster)]
  plot = plot + scale_colour_manual(values = c('grey90', cluster_color))
  #if(cortical){plot = plot + 
  #  geom_hline(yintercept=boundary_dorsal, linetype = "dashed",color = cluster_color) +
  #  geom_hline(yintercept=boundary_ventral, linetype = "dashed",color = cluster_color)}
  return(plot)
}

plot_clusters_vertical_spatial_no_grid_merge <- function(sobj, stream, clustering = NULL, anterior = FALSE, 
                                                   cluster_color =  '#CB2A55', pt.size = 1, space = "H", 
                                                   arc = TRUE, cortical = TRUE, x_width = 35, force_idents = NA)
{
  if(stream == 'CC'){
    ixs = grepl("CC", sobj$area)
  } else if (stream == 'PirC'){
    ixs = grepl("TC", sobj$area) & grepl("164", sobj$area)
  } else if (stream == 'TC'){
    ixs = grepl("TC", sobj$area) & grepl("408", sobj$area)
  } else{
    stop("Invalid stream name provided")
  }
  
  sobj = sobj[, ixs]
  final_df <- data.frame()
  x_adj = 5
  if(any(is.na(force_idents))){
    idents = levels(Idents(sobj))
  } else{
    idents = force_idents
  }
  #x_width = ifelse(get_slice(sobj) == '164', )
  for(cluster in idents){
    cluster_identity = Idents(sobj) == cluster
    coordinates <- Embeddings(sobj, reduction = space)
    gene_df <- as.data.frame(cbind(coordinates, cluster_identity))
    colnames(gene_df) <- c('X', 'Y', 'clust')
    if(arc){
      dms_ixs = which(sobj$image %in% dms | sobj$image %in% cc)
      gene_df$Y[dms_ixs] <- gene_df$Y[dms_ixs] + 50
    }
    
    if(cortical){
      if(get_slice(sobj) == '164'){
        tc = antr_TC_Cx
        cc = antr_CC_Cx
      } else{
        tc = post_TC_Cx
        cc = post_CC_Cx
      }
      non_tc_ix = which(!sobj$image %in% tc)
      cc_ix = which(sobj$image %in% cc)
      gene_df$Y[non_tc_ix] <- gene_df$Y[non_tc_ix] + 15
      gene_df$Y[cc_ix] <- gene_df$Y[cc_ix] + 15
      gene_df$cluster = cluster
      #min_vms <- min(gene_df$Y[which(sobj$image %in% vms)])
      #max_tc <- max(gene_df$Y[which(sobj$image %in% tc)])
      #max_dms <- max(gene_df$Y[which(sobj$image %in% dms)])
      #min_cc <- min(gene_df$Y[which(sobj$image %in% cc)])
      #boundary_dorsal <- mean(min_vms, max_tc)
      #boundary_ventral <- mean(max_dms, min_cc)
    }
    
    gene_df$X = gene_df$X + x_adj
    x_adj = x_adj + x_width
    final_df <- rbind(final_df, gene_df)
  }
  
  ix_164 = startsWith(rownames(final_df), '164')
  ix_408 = startsWith(rownames(final_df), '408')
  
  min_164 = min(final_df$Y[ix_164])
  min_408 = min(final_df$Y[ix_408])
  
  max_164 = max(final_df$Y[ix_164])
  max_408 = max(final_df$Y[ix_408])
  max_avg = (max_164 + max_408) / 3
  
  final_df$Y[ix_164] = final_df$Y[ix_164] - min_164
  final_df$Y[ix_408] = final_df$Y[ix_408] - min_408
  
  final_df$Y[ix_164] = (final_df$Y[ix_164] / (max_164 - min_164)) * max_avg
  final_df$Y[ix_408] = (final_df$Y[ix_408] / (max_408 - min_408)) * max_avg
  
  #final_df$Y[ix_164] = 0
  
  final_df <- final_df %>% dplyr::arrange(clust)
  
  final_df$cluster[which(final_df$clust == 0)] = NA
  
  plot <- ggplot(final_df, aes(x = X, y = Y, color = factor(cluster, levels = idents))) + geom_point(size = pt.size, alpha = 1) +  
    theme_classic() + #ggtitle(cluster) + 
    NoAxes() + NoLegend() + 
    coord_fixed(ratio = 0.5)  + theme(title = element_text(face = 'bold', size = rel(0.8), hjust = 1)) 
  #cluster_color = scales::hue_pal()(length(idents))[which(idents == cluster)]
  if(length(levels(sobj)) != length(idents)){
    idents = idents[which(idents %in% levels(sobj))]
  }
  plot = plot + scale_colour_manual(values = get_cluster_colors(idents), na.value = 'grey90')
  #if(cortical){plot = plot + 
  #  geom_hline(yintercept=boundary_dorsal, linetype = "dashed",color = cluster_color) +
  #  geom_hline(yintercept=boundary_ventral, linetype = "dashed",color = cluster_color)}
  return(plot)
}

