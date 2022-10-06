library(tidyverse)
library(patchwork)
library(DEqMS)
library(fgsea)

##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/')
protein_flattened_wide = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed.tsv')
protein_long_meta = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed_long_w_meta.tsv')
phos_flattened_wide = read_tsv('phospho/processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed.tsv')
phos_long_meta = read_tsv('phospho/processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed_long_w_meta.tsv')
meta = read_tsv('../Data/Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')
markers = read_tsv('../Resource/Marker/KEGG_PI3K-AKI_genelist.txt')
markers_pik3ca = read_tsv('../Resource/Marker/PI3K_score.tsv')
# https://www.genome.jp/pathway/hsa04151


##############################
## Source functions
walk(list.files('script/proteomics/src/', full.names = T), source)

####################################################################################################
### Clean data
####################################################################################################
## Marker
markers_pik3ca_clean = markers_pik3ca %>% mutate(Phosphosite = str_remove(Phosphosite, '^p') %>% str_remove('_')) %>% 
  select(Target_type, Gene, Phosphosite, `Type of regulation`, Pathway)

## Phos
phos_cols =c('Gene', 'LocalizedPhosphosites', 'New SampleID', 'Phospho_Site_level',# data
  'Histology', 'Passage (pro)', 'Treatment', 'PIK3CA mutation') # metadata
phos_long_clean = phos_long_meta %>% filter(Species == 'HUMAN') %>% 
  select(all_of(phos_cols)) %>% dplyr::rename(Phosphosite = LocalizedPhosphosites)

## Pro
pro_cols = c('Protein', 'New SampleID', 'Protein_level',# data
             'Histology', 'Passage (pro)', 'Treatment', 'PIK3CA mutation') # metadata
pro_long_clean = protein_long_meta %>% filter(Species == 'HUMAN') %>% 
  select(all_of(pro_cols)) 

# 
# ####################################################################################################
# ### Plot PIK3CA Phospho score
# ####################################################################################################
# phos_marker_pik3ca = markers_pik3ca_clean %>% filter(Target_type =='Phospho')
# phos_pik3ca_level =  inner_join(phos_marker_pik3ca, phos_long_clean, by = c('Gene','Phosphosite')) 
# 
# # Check all possible phososite for genes
# all_positive_df = phos_long_clean %>% filter(Gene %in% phos_marker_pik3ca$Gene) %>% select(Gene, Phosphosite) %>% distinct 
# 
# phos_long_clean   %>% filter(Gene == 'AKT1',Phosphosite == 'T308') %>% dim
# phos_pik3ca_level %>% filter(Gene == 'AKT1',Phosphosite == 'T308') %>% dim
# 
# 
# 
# 
# asdfasdf





####################################################################################################
### Plot PIK3CA PRO score
####################################################################################################
pro_marker_pik3ca_df = markers_pik3ca_clean %>% filter(Target_type =='Protein')
pro_marker_pik3ca_gene = pro_marker_pik3ca_df$Gene
####################################################################################################
### Plot top variance
## turn into function? For easier Pan-can and Inidiviaul can plotting
####################################################################################################
protein_mtx_human = protein_flattened_wide %>% filter(Species == 'HUMAN') %>% 
  select(Protein, contains('PA0')) %>% as.data.frame %>% 
  column_to_rownames('Protein') 

# Marker only
pathways = c("Hormone B Signaling pathway", "PI3K/Akt Signaling pathway", 
             "Apoptosis Signaling pathway", "DNA Damage Signaling pathway", 
             "Cell Cycle Signaling pathway", "Hormone A Signaling pathway")
selected_genes = pro_marker_pik3ca_df %>% filter(Pathway %in%  pathways[4:5]) %>% pull(Gene)
protein_mtx_marker_human = protein_mtx_human %>% .[intersect(rownames(.), selected_genes),]# Marker only

########
# Palette 
#########
col_anno = c("Institution", "Passage (pro)", "Treatment", 'Histology', "PIK3CA mutation")
palette_list_use = MakeMetaPalette(meta = meta, columns_use = col_anno, existing_palette_list = palette_meta)
palette_list_pancan = c(palette_list_use)



#####################################################
#####################################################
#### PANCAN
# Step1: Make Plot object 
pancan_obj = MakeDataObj(protein_mtx_marker_human, meta = meta, top_var = 0, na_cutoff = 0.5, filter_var = F) %>%  # set top_var = 0 to keep all genes
    MakeCluster(k=3) %>% 
    MakeUmap()
# Step 2: Make heatmap
# heatmap
dir.create('proteomic/figure/Pancan/PIK3CA_Pathway/')
pdf('proteomic/figure/Pancan/PIK3CA_Pathway/Heatmap.pdf', width = 10, height = 6)
MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column = F, show_row_names = T,
                  anno_pch_cols = c('Histology', "Treatment"),
                  anno_regular = c("Institution","PIK3CA mutation","Passage (pro)", 'cluster'))
MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column_by = 'Histology', show_row_names = T)
MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column_by = 'cluster',
                  anno_pch_cols = c('Histology', "Treatment"),
                  anno_regular = c("Institution","PIK3CA mutation","Passage (pro)", 'cluster'))

MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column_by = 'PIK3CA mutation',
                  anno_pch_cols = c('Histology', "Treatment"),
                  anno_regular = c("Institution","PIK3CA mutation","Passage (pro)", 'cluster'))
dev.off()
pdf('proteomic/figure/Pancan/PIK3CA_Pathway/Umap.pdf', width = 10, height = 10)
MakeMultipleDimPlot(pancan_obj, palette_list = palette_list_use)
dev.off()

Heatmap(column_title_gp = gpar(col = 'red'))
#####################################################
### INdividual
## RCC 
obj_list[['RCC']] = MakeDataObj(protein_mtx_human, meta, histology = 'RCC') %>% MakeCluster(k=4) %>% MakeUmap
dir.create('figure/RCC/')
pdf('figure/RCC/Heatmap.pdf', width = 10, height = 10)
MakeTopProHeatmap(obj_list[['RCC']], palette_list = palette_list_pancan, split_column = F)
MakeTopProHeatmap(obj_list[['RCC']], palette_list = palette_list_pancan, split_column_by = 'PIK3CA mutation')
MakeTopProHeatmap(obj_list[['RCC']], palette_list = palette_list_pancan, split_column_by = 'cluster')
dev.off()
pdf('figure/RCC/Umap.pdf', width = 10, height = 10)
MakeMultipleDimPlot(obj_list[['RCC']], palette_list = palette_list_pancan, column_names = c('PIK3CA mutation','Treatment','cluster'))  
dev.off()

# Breast 
obj_list[['Breast']] = MakeDataObj(protein_mtx_human, meta, histology = 'breast cancer') %>% MakeCluster(k=4) %>% MakeUmap
dir.create('figure/Breast/')
pdf('figure/Breast/Heatmap.pdf', width = 10, height = 10)
MakeTopProHeatmap(obj_list[['Breast']], palette_list = palette_list_pancan, split_column = F) 
MakeTopProHeatmap(obj_list[['Breast']], palette_list = palette_list_pancan, split_column_by = 'PIK3CA mutation') 
MakeTopProHeatmap(obj_list[['Breast']], palette_list = palette_list_pancan, split_column_by = 'cluster')
dev.off()
pdf('figure/Breast/Umap.pdf', width = 7, height = 7)
MakeMultipleDimPlot(obj_list[['Breast']], palette_list = palette_list_pancan, column_names = c('PIK3CA mutation','Treatment','cluster'))  
dev.off()


## Others
histology_plt = c("CRC", "Lung Adenocarcinoma", "melanoma")
for(hist in histology_plt){
  obj_list[[hist]] = MakeDataObj(protein_mtx_human, meta, histology = hist) %>% MakeCluster(k=2) %>% MakeUmap
  dir.create(str_glue('figure/{hist}/'))
  pdf(str_glue('figure/{hist}/Heatmap.pdf'), width = 10, height = 10)
  MakeTopProHeatmap(obj_list[[hist]], palette_list = palette_list_pancan, split_column = F) %>% draw
  MakeTopProHeatmap(obj_list[[hist]], palette_list = palette_list_pancan, split_column_by = 'PIK3CA mutation') %>% draw
  MakeTopProHeatmap(obj_list[[hist]], palette_list = palette_list_pancan, split_column_by = 'cluster') %>% draw
  dev.off()
  pdf(str_glue('figure/{hist}/Umap.pdf'), width = 7, height = 7)
  MakeMultipleDimPlot(obj_list[[hist]], palette_list = palette_list_pancan, column_names = c('PIK3CA mutation','Treatment','cluster'))%>% plot
  dev.off()
}




########################################
## PICK3CA activation score
#######################################
pathway_table = read_tsv('../../Resource/pathway/PI3K_pathway_protein_phoso.tsv')

pathway_table %>% mutate(Exist = ifelse(Protein %in% protein_long_meta$Protein, 'Exist','No')) %>% 
  count(Pathway, Exist) %>% 
  pivot_wider(id_cols = Pathway, names_from = 'Exist', values_from = 'n')


#### Consensus?

