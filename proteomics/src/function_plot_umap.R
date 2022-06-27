

################ 
### Plot
################
## Plot Umap result
MakeDimPlot = function(data_obj, group.by = 'cluster'){
  if(data_obj@coordinate %>% length ==0) stop('Missing UMAP coordinate. Run MakeUmap first')
  umap_df = data_obj@coordinate$umap_table
  
  plot_table = bind_cols(umap_df, data_obj@meta.data[rownames(umap_df),])
  plot_table %>% ggplot(aes(x = x ,y =y, color = .data[[group.by]] )) + 
    geom_point() + 
    labs(title = group.by) +
    cowplot::theme_cowplot()+ 
    theme(aspect.ratio = 1)
}


# Dim
MakeMultipleDimPlot = function(data_obj, column_names = c('Histology','PIK3CA mutation', 'cluster'), palette_list){
  p_umap_list = map(column_names, ~MakeDimPlot(data_obj, group.by =.)) %>% setNames(column_names)
  # palette
  if(!missing(palette_list)){
    column_exist = intersect(column_names, names(palette_list))
    for(column in column_exist){
      palette_use = palette_list[[column]][names(palette_list[[column]]) %in% data_obj@meta.data[[column]]]
      p_umap_list[[column]] =  p_umap_list[[column]] + scale_color_manual(values = palette_use)
    }
  }
  wrap_plots(p_umap_list, guides = 'collect') & 
    theme(
      # Boarder
      panel.background = element_rect(color = 'gray30', size = 2),
      # Title
      plot.title = element_text(hjust = 0.5), 
          # remove axis:
          axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), 
          # legends
          legend.position = 'bottom', 
          legend.box = "vertical",
          legend.box.just  = 'left')
}

