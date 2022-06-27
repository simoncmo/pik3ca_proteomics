library(ComplexHeatmap)
#########################
## Heatmap
####################
# Palette 
# Idea for helper functionUse 'multiple palette to generate color palettes - 
# A. from this https://github.com/EmilHvitfeldt/r-color-palettes 
# B. Use RColorBrewer
rcolorqual = function(n=74){ 
  # All qual -  max 74 colors
  qual_pals_table = RColorBrewer::brewer.pal.info %>% filter(category =='qual') %>% rownames_to_column('pal_name')
  # Get all col
  all_pal_list = pmap(qual_pals_table, function(pal_name, maxcolors, ...) RColorBrewer::brewer.pal(n = maxcolors, name = pal_name)) #%>% setNames(qual_pals_table %>% rownames)
  return(unlist(all_pal_list)[1:n]) 
}

## Palette
make_palette = function(vector){
  items = unique(vector)
  rcolorqual() %>% sample(length(items)) %>% setNames(items)
}


## Generate Palette for Meta
MakeMetaPalette = function(meta, columns_use, existing_palette_list = list()){
  if(missing(columns_use)) columns_use = names(meta) # All columns
  no_palette_col = setdiff(columns_use, names(existing_palette_list)) # Missing palette
  meta = as.data.frame(meta)
  if(length(no_palette_col) ==1) col_item_list = list(unique(meta[, no_palette_col] )) %>% setNames(no_palette_col)
  else col_item_list = apply(meta[, no_palette_col, drop=F], 2, unique) # Get all items per column
  #return(col_item_list)
  # Assign palette and add to list
  new_palette_list = c(map(col_item_list, make_palette), existing_palette_list)
  # remove NA
  new_palette_list = map(new_palette_list, ~.[!is.na(names(.))])
  return(new_palette_list)
}

### Legend for pch too annot
MakeLegendList = function(palette_list, groups_use){
  palette_list = if(!missing(groups_use)) palette_list[groups_use] else palette_list
  imap(palette_list, function(pal, name){
    lgd = Legend(labels = pal %>% names, legend_gp = gpar(fill = pal), title = name)
  })
}

### Top annotation
MakeTopAnnotation = function(anno_pch_cols, anno_regular, meta_data, palette_list){
  # Palette
  col_column_select = c(anno_pch_cols, anno_regular)
  col_palette_list = MakeMetaPalette(meta = meta_data, columns_use = col_column_select, existing_palette_list = palette_list)
  
  ###### 
  # TOP ANNOTATIONs
  # Create anno_simple
  # pch
  anno_pch_list = map2(meta_data[,anno_pch_cols, drop=F], col_palette_list[anno_pch_cols], function(data, pal){
    anno_simple(data, pch =  data, col = pal)
  })

  # regular
  anno_regular_list = meta_data[,anno_regular] %>% as.list()
  anno_regular_col_list = imap(col_palette_list[anno_regular], function(pal, column){
    column = pal
  })
  anno_regular_list = c(anno_regular_list, list(col = anno_regular_col_list))
  # combine
  anno_arg_list = c( anno_regular_list, anno_pch_list)
  
  # Create top anno
  top_annot = do.call(HeatmapAnnotation, anno_arg_list)
  return(top_annot)
}


### Heatmap func
MakeTopProHeatmap = function(data_obj, na_cutoff = 0.5, top_var = 0.9, histology, max_gene_plt = 3000,
                             anno_pch_cols = c( "Treatment"),
                             anno_regular = c("Institution",'Histology', "PIK3CA mutation","Passage (pro)", 'cluster'),
                             palette_list = list(),
                             split_column = T,
                             split_column_by = 'cluster',
                             gene_highlight = markers$Gene, # By default highlight PIK3CA genes
                             ...
){
  
  # Select data
  top_genes = data_obj@gene.meta %>% filter(Top_gene == TRUE) %>% rownames
  #return(data_obj@data)
  mtx_selected = data_obj@data[top_genes, ]
  meta_data = data_obj@meta.data
  
  # # Check if filter
  # if(nrow(mtx_selected) > max_gene_plt){ # Some wheird about ths RANK ... PAUSED
  #   message(str_glue('Totally {nrow(mtx_selected)} genes. Plot top {max_gene_plt} instead. Change max_gene_plt if needed'))
  #   top_n_genes  = data_obj@gene.meta %>% filter(Var_rank <max_gene_plt) %>% rownames 
  #   mtx_selected = mtx_selected[top_n_genes, ]
  # }
  ############
  # Visualize - Top 10% Most variable genes 
  ############
  # Palette
  col_column_select = c(anno_pch_cols, anno_regular)
  col_palette_list = MakeMetaPalette(meta = meta_data, columns_use = col_column_select, existing_palette_list = palette_list)

  # ###### 
  # # TOP ANNOTATIONs
  top_annot = MakeTopAnnotation(anno_pch_cols = anno_pch_cols, anno_regular = anno_regular, meta_data = meta_data, palette_list = palette_list)
  # Top anno legend for pch
  anno_pch_lgd_list = MakeLegendList(palette_list[anno_pch_cols])
  
  #### ROW ANNOTATION
  # Gene Mark annotation
  gene_label = intersect(markers$Gene, rownames(mtx_selected))
  rw_anno = if(length(gene_label)>0){
    rowAnnotation(foo = anno_mark(at = match(gene_label, rownames(mtx_selected)), labels = gene_label))
  }else{
    rowAnnotation(foo = anno_empty())
  }
  
  #return(mtx_selected)
  heatmap_param = c(list(matrix = mtx_selected %>% t %>% scale %>% t,
                       right_annotation = rw_anno, 
                       top_annotation = top_annot, 
                       name = 'Protein\nLevel\nScaled',
                       show_row_names = F,
                       show_row_dend = F,
                       show_column_names = F), list(...))
  
  # split
  if(split_column) heatmap_param = c(heatmap_param, list(column_split = meta_data[[split_column_by]]))
  hm = do.call(Heatmap, args = heatmap_param)
  
  # Add legends from pch
  hm = draw(hm, annotation_legend_list = anno_pch_lgd_list)

}
