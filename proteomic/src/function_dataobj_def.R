library(tidyverse)
####################
## Define data container
####################
# S4 define class
dataobj <- setClass('dataobj', slots = list(data = 'matrix', meta.data = 'data.frame', gene.meta = 'data.frame', settings = 'list', cluster = 'list', coordinate = 'list'))
# Set show function
setMethod('show','dataobj', function(object){
  cat('Size of data:', toString(dim(object@data)),'\n')
  cat('Number of sample:', nrow(object@meta.data),'\n')
  cat('Meta columns:', toString(names(object@meta.data)),'\n')
  cat('Setting:', toString(unlist(object@settings)),'\n')
}) 
## Make data frame 
MakeDataObj = function(mtx, meta, na_cutoff = 0.5, top_var = 0.9, histology, filter_na = T, filter_var = T
){
  
  #### Meta and Sample filtering
  meta_data = meta %>% filter(!str_detect(`Sample aliquot ID`,'reference')) %>% as.data.frame
  rownames(meta_data) = meta_data$`Sample aliquot ID`
  if(!missing(histology)) meta_data = meta_data %>% filter(Histology == histology)  #c(NA, "RCC", "CRC", "Lung Adenocarcinoma", "melanoma", "breast cancer")
  sample_use = meta_data$`Sample aliquot ID`
  mtx = mtx[,sample_use]
  
  #### GENE filtering
  # Filter by NA percent
  if(filter_na){
    pro_na_percent = mtx %>% apply(1, is.na) %>% colMeans()
    pro_filtered   = pro_na_percent[pro_na_percent < na_cutoff] %>% names
    mtx_filtered   = mtx[pro_filtered, ]
  }
  
  # Selecy by var
  pro_var = mtx_filtered %>% apply(1, var, na.rm=T)
  var_cutoff = quantile(pro_var, top_var) # top 10 most variable genes 
  pro_select = if(filter_var) pro_var[pro_var > var_cutoff] %>% names else names(pro_var)
  mtx_selected = mtx_filtered[pro_select,]
  
  # meta for gene
  gene.meta = data.frame(Var = pro_var, Var_rank = rank(pro_var), Top_gene = pro_var > var_cutoff)
  
  
  # Make obj
  settings = list(na_cutoff = na_cutoff, top_var=top_var)
  obj = new('dataobj', 
            # Need to decide how to strore original and FILTERED mtx
            data = mtx_selected %>% as.matrix, #mtx_filtered %>% as.matrix, 
            gene.meta = gene.meta,
            meta.data = meta_data,
            settings = settings)
  return(obj)
}



#########################
## Make cluster and Umap
## Use list to enclose data obj
####################
library(umap)
### Caculate cluster - kmeans
MakeCluster =function(data_obj, k = 5){
  cat('Using kmeans with k = ',k ,'\n')
  kmean_data = data_obj@data %>% t; kmean_data[is.na(kmean_data)] =0 # Need replace na with 0
  kmeans_result = kmeans(kmean_data ,centers = k)
  data_obj@meta.data = data_obj@meta.data %>% mutate(cluster = kmeans_result$cluster[rownames(.)] %>% as.character()) # Add to meta
  data_obj@cluster = list(kmeans = kmeans_result)
  return(data_obj)
}
### Caculate umap
MakeUmap = function(data_obj, seed = 42){
  set.seed(seed)
  umap_data = data_obj@data %>% t ; umap_data[is.na(umap_data)] =0 # Need replace na with 0
  # set n_neighbor
  n_neighbors = if (nrow(umap_data) < 15) nrow(umap_data)-1 else 15
  if(n_neighbors<=1){ # too few sample to cluster
    warnings('Less and equal to 2 sample. Too few to cluster')
    return(data_obj)
  }
  umap_out = umap::umap(d = umap_data, n_neighbors = n_neighbors)
  
  umap_table = umap_out$layout %>% as.data.frame %>% setNames(c('x','y'))
  data_obj@coordinate = list(umap_out = umap_out, umap_table = umap_table)
  return(data_obj)
}
