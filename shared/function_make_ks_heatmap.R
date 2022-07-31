##############################
# Created 7/30/2022
# Simon Mo
# For PIK3CA PDX 40 sample Only at the moment 
##############################
# Main Heatmap

# DEMO:
# MakeKSHeatmap(protein_flattened_wide, phos_flattened_wide, meta_use, kinase_table = kinase_table, 
#               kinase_of_interest = 'EGFR',
#               na_cutoff_kinase = 1, na_cutoff_substrate = 0.75)
##############################

MakeKSHeatmap = function(protein_wide_mtx, phos_wide_mtx, meta_table, kinase_table, kinase_of_interest=NA,
                        na_cutoff_kinase = 1, na_cutoff_substrate = 0.4,
                            # Palette
                        anno_pch_cols = c("Treatment"),
                        anno_regular = c("Institution",'Histology', "PIK3CA mutation","Passage (pro)"),
                        palette_list = palette_list_pancan
                    ){
    hm_kinase = MakeKSHeatmap_base(protein_wide_mtx, phos_wide_mtx, meta_table, kinase_table, kinase_of_interest, mode ='Kinase',
                        na_cutoff = na_cutoff_kinase, 
                        sample_order = NULL,
                        anno_pch_cols, anno_regular, palette_list)
    # note: substrate sample order depends on kinase
    kinase_sample_order = hm_kinase@matrix %>% colnames
    hm_substrate = MakeKSHeatmap_base(protein_wide_mtx, phos_wide_mtx, meta_table, kinase_table, kinase_of_interest, mode ='Substrate',
                        na_cutoff = na_cutoff_substrate,
                        sample_order = kinase_sample_order,
                        anno_pch_cols, anno_regular, palette_list)
    hm_kinase %v% hm_substrate

}

MakeKSHeatmap_base = function(protein_wide_mtx, phos_wide_mtx, meta_table ,kinase_table, kinase_of_interest=NA, mode = c('Kinase','Substrate'),
                        na_cutoff,
                        # sample order 
                        sample_order = NULL,
                        # Palette
                        anno_pch_cols, anno_regular, palette_list
                    ){
    # Parameters
    if(!kinase_of_interest %in% protein_wide_mtx$Protein) stop(str_glue("{kinase_of_interest} not found in the Protein table provided. Please Double check"))
    mode = match.arg(mode)
    if(mode == 'Substrate') substrate_vec = kinase_table %>% filter(Kinase == kinase_of_interest) %>% pull(Substrate) %>% setdiff(., kinase_of_interest)
    gene_use = if(mode == 'Substrate') substrate_vec else kinase_of_interest

    # MATRIX: Construct
    pro_table = protein_wide_mtx %>% 
        filter(Species == 'HUMAN') %>% 
        filter(Protein %in% all_of(gene_use)) %>% 
        mutate(Feature = Protein)

    # Using Idex to rename to make sure not duplications
    phos_table = phos_wide_mtx %>% 
      filter(Species == 'HUMAN') %>% 
      filter(Gene %in% all_of(gene_use)) %>% 
      MakeDistinctPhosFeature()

    # MATRIX: Combine Pro and Phos matrix
    combined_matrix = bind_rows(
        pro_table %>% select(Feature, all_of(sample_ids)),
        phos_table %>% select(Feature, all_of(sample_ids))
    ) %>% column_to_rownames('Feature')

    # MATRIX:ROW: Filter GENE by NA ratio
    # Remove gene "ABOVE" the cutoff
    feature_keep = rownames(combined_matrix)[rowMeans(is.na(combined_matrix)) <= na_cutoff]
    combined_matrix = combined_matrix[feature_keep, ]

    # MATRIX: Fill na - with 0 
    combined_matrix[is.na(combined_matrix)] = 0

    # MATRIX:COLUMN:ORDER Get order of the "Kinase protein"
    # For kinase only
    if(is.null(sample_order)){
        if(mode == 'Kinase'){
        sample_order = combined_matrix %>% .[kinase_of_interest,] %>% 
            #select(sample_ids) %>% 
            t %>% as.data.frame %>% setNames('Protein_level')  %>%
            arrange(desc(Protein_level)) %>% 
            rownames()        
        }else{
            stop('In Substrate mode, please provide sample order')
        }
    }else{
        message('Sample order from Kinase')
    }

    # MATRIX:ROW:ORDER Get sample order without Kinase protein
    if(mode == 'Kinase'){
        # Kinase only: Protein on TOP
        hclust_result = combined_matrix %>% .[rownames(.)!=kinase_of_interest, ] %>% dist %>% hclust # %>% .$order
        feature_ordered = hclust_result$labels[hclust_result$order] %>% rev
        feature_ordered = c(kinase_of_interest, feature_ordered) # Add kinase protein back
    }else{
        hclust_result = combined_matrix %>% dist %>% hclust # %>% .$order
        feature_ordered = hclust_result$labels[hclust_result$order] %>% rev
    }
    
    
    # MATRIX:COLUMN:ORDER: Reorder matrix and meta.data
    # SAMPLE_ANNOT:ORDER
    combined_matrix = combined_matrix[,sample_order]
    meta_use = meta_table[sample_order,]


    # MATRIX:ROW: Color PROTEIN with 'red' 
    # MATRIX:COL: Color Sample by Cancer Type
    row_text_color_vec = GetGeneTextColor(combined_matrix)[feature_ordered]
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

    # Make heatmap
    p_heatmap = Heatmap(combined_matrix[feature_ordered, sample_order] %>% as.matrix, 
            cluster_columns = F, 
            cluster_rows = F,
            top_annotation = if(mode == 'Kinase') top_annot else NULL, 
            row_names_gp = gpar(col = row_text_color_vec),
            column_names_rot = -90,
            column_names_gp = gpar(fill = col_text_color_vec),
            right_annotation = MakeKSRowBlock(mode), 
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


