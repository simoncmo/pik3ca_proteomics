##############################
# Created 8/26/2022
# Simon Mo
# For PIK3CA PDX 40 sample Only at the moment 
##############################
# Main Heatmap

# DEMO:
# MakeKSHeatmap(protein_flattened_wide, phos_flattened_wide, meta_use, kinase_table = kinase_table, 
#               kinase_of_interest = 'EGFR',
#               na_cutoff_kinase = 1, na_cutoff_substrate = 0.75)
##############################

MakeProPhosHeatmap = function(protein_wide_mtx, phos_wide_mtx, gene_use, meta_table,
                        na_cutoff = 0.5,
                        sort_row = T, sort_col =T,
                        split_by_gene = F,
                        # Palette
                        anno_pch_cols = c("Treatment"),
                        anno_regular = c("Institution",'Histology', "PIK3CA mutation","Passage (pro)"),
                        palette_list =palette_list_pancan
                    ){
    # MATRIX: Combine Pro and Phos matrix
    combined_matrix = MakeCombinedProPhosMtx(protein_wide_mtx, phos_wide_mtx, gene_use)

    # MATRIX:ROW: Filter GENE by NA ratio and FillNA
    combined_matrix = FilterMtxAndFillNa(combined_matrix, na_cutoff)

    # MATRIX:COLUMN:ORDER 
    # MATRIX:ROW:ORDER 
    combined_matrix = SortMtxHclust(combined_matrix, sort_row = sort_row, sort_col = sort_col)
    
    # MATRIX:COLUMN:ORDER: Reorder meta.data
    # SAMPLE_ANNOT:ORDER
    meta_use = meta_table[colnames(combined_matrix),]

    # MATRIX:ROW: Color PROTEIN with 'red' 
    # MATRIX:COL: Color Sample by Cancer Type
    row_text_color_vec = GetGeneTextColor(combined_matrix)
    col_text_color_vec = GetSampleTextColor(meta_use, palette_list, 'Histology')
    

    ############
    # TOP ANNOTATIONs: Parameter
    ############
    # Palette
    col_column_select = c(anno_pch_cols, anno_regular)
    col_palette_list = MakeMetaPalette(meta = meta_use, columns_use = col_column_select, existing_palette_list = palette_list)
    # TOP ANNOTATIONs: CREATE
    top_annot = MakeTopAnnotation(anno_pch_cols = anno_pch_cols, anno_regular = anno_regular, meta_data = meta_use, palette_list = palette_list)
    # Top anno legend for pch
    # anno_pch_lgd_list = MakeLegendList(palette_list[anno_pch_cols])


    # MATRIX: Color
    col_fun = colorRamp2(seq(-4,4, length.out=11), RColorBrewer::brewer.pal(n=11, name = 'RdBu'))
    col_fun(seq(-3, 3))

    # ANNO: RIGHT
    m = as.matrix(combined_matrix)
    anno_right = rowAnnotation(foo = anno_boxplot(m, height = unit(1, "cm"), outline = FALSE))

    # Make heatmap
    p_heatmap = Heatmap(m, 
            cluster_columns = F, 
            cluster_rows = F,
            row_split = if(!split_by_gene) NULL else rownames(combined_matrix) %>% str_remove('_.+'),
            top_annotation = top_annot, 
            right_annotation = anno_right,
            row_names_gp = gpar(col = row_text_color_vec),
            column_names_rot = -90,
            column_names_gp = gpar(fill = col_text_color_vec),
            name = 'Protein\nPhosphoprotein\nLevel',
            col = col_fun) 
    
    p_heatmap
}




##############################
# UTILITY FUNCTIONS
##############################
# [Function]: Make Label block on the right 
MakeKSRowBlock = function(group = c('Kinase','Substrate')){
  group = match.arg(group) 
  
  # Block Color
  col_block = RColorBrewer::brewer.pal(n=3, 'Dark2')[1:2] %>% setNames(c('Substrate','Kinase'))
  col_block = col_block[group]
  
  rowAnnotation(foo = anno_block(labels_gp = gpar(col = 'white'),
                                 gp = gpar(fill = col_block, col = 'transparent'),
                                 labels_rot = -90,
                                 labels = group,
                                 width = unit(15, units = 'pt')))
}


# [Function]: Matrix:ROW:Text:Color
GetGeneTextColor = function(mtx, col_phos_pro){
  if(missing(col_phos_pro)) col_phos_pro = RColorBrewer::brewer.pal(7, 'Set1')[3:4] %>% setNames(c('Phos','Pro'))
  row_text_color_vec = ifelse(mtx %>% rownames %>% str_detect('_'), col_phos_pro[[1]], col_phos_pro[[2]])
  row_text_color_vec %>% setNames(mtx %>% rownames)
}

# [Function]: Matrix:ROW:Text:Color
GetSampleTextColor = function(meta_table, palette_list, column_use){
  if(missing(column_use)){
    column_use = palette_list %>% names %>% .[1]
    message(str_glue('Missing column_use, use first item in the palette_list: {column_use}'))}
  if(missing(meta_table)| missing(palette_list)) stop('Missing meta_table or palette')
  # Choose column and palette to use 
  palette_list[[column_use]][meta_table[[column_use]]]
}


# Utility function: Make Unique Phos name
MakeDistinctPhosFeature = function(phos_wide_mtx){
    # Using Idex to rename to make sure not duplications
    phos_table = phos_wide_mtx %>% 
      mutate(Feature = str_c(Gene, FirstPossiblePhosphoSite, SecondPossiblePhosphoSite, LocalizedPhosphosites, sep='_')) %>% 
      mutate(Feature_short = str_c(Gene, LocalizedPhosphosites, sep='_'))

    # MATRIX:ROW: Replace Feature with 'short' version if distinct
    distinct_phospho_name = phos_table %>% count(Feature_short) %>% filter(n == 1) %>% pull(Feature_short)
    phos_table %>% mutate(Feature = ifelse(Feature_short %in% distinct_phospho_name, 
                                                                 Feature_short, 
                                                                 Feature))
}

##############
# Matrix combine and utility function

# Step1 : Combine
MakeCombinedProPhosMtx = function(protein_wide_mtx, phos_wide_mtx, gene_use){

  sample_ids = protein_wide_mtx %>% names %>%  str_subset('PA')
    # MATRIX: Construct
    pro_table = protein_wide_mtx %>% 
        filter(Species == 'HUMAN') %>% 
        mutate(Feature = Protein)

    # Using Idex to rename to make sure not duplications
    phos_table = phos_wide_mtx %>% 
      filter(Species == 'HUMAN') %>%
      MakeDistinctPhosFeature()

    # Select gene
    if(!missing(gene_use)){
      pro_table = pro_table %>% filter(Protein %in% all_of(gene_use)) 
      phos_table = phos_table %>% filter(Gene %in% all_of(gene_use)) 
    }

    # MATRIX: Combine Pro and Phos matrix
    combined_matrix = bind_rows(
        pro_table %>% select(Feature, all_of(sample_ids)),
        phos_table %>% select(Feature, all_of(sample_ids))
    ) %>% column_to_rownames('Feature')
    return(combined_matrix)
}

# Step2 : Process NA
FilterMtxAndFillNa = function(combined_matrix, na_cutoff = 0.5){
  # MATRIX:ROW: Filter GENE by NA ratio
    # Remove gene "ABOVE" the cutoff
    feature_keep = rownames(combined_matrix)[rowMeans(is.na(combined_matrix)) <= na_cutoff]
    combined_matrix = combined_matrix[feature_keep, ]

    # MATRIX: Fill na - with 0 
    combined_matrix[is.na(combined_matrix)] = 0
    return(combined_matrix)
}

# Step3 : sort with hclust
SortMtxHclust = function(mtx, sort_row = T, sort_col = T){
  # Row
  if(sort_row){
    row_order = mtx %>% dist %>% hclust %>% {.$labels[.$order]}
    mtx = mtx[row_order, ]
  }
  if(sort_col){
    col_order = mtx %>% t %>% dist %>% hclust %>% {.$labels[.$order]}
    mtx = mtx[, col_order]
  }
  
  return(mtx)
}