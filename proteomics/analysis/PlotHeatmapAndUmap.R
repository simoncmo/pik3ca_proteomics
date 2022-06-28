library(tidyverse)
library(patchwork)
library(DEqMS)
library(fgsea)

##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/')
protein_flattened_wide = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed.tsv')
protein_long_meta = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed_long_w_meta.tsv')
meta = read_tsv('../Data/Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')
markers = read_tsv('../Resource/Marker/KEGG_PI3K-AKI_genelist.txt')
# https://www.genome.jp/pathway/hsa04151


##############################
## Source functions
walk(list.files('script/proteomics/src/', full.names = T), source)


####################################################################################################
### Plot top variance
## turn into function? For easier Pan-can and Inidiviaul can plotting
####################################################################################################
protein_mtx_human = protein_flattened_wide %>% filter(Species == 'HUMAN') %>% 
  select(Protein, contains('PA0')) %>% as.data.frame %>% 
  column_to_rownames('Protein') 

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
pancan_obj = MakeDataObj(protein_mtx_human, meta = meta, top_var = 0.8) %>% 
    MakeCluster(k=5) %>% 
    MakeUmap()
# Step 2: Make heatmap
# heatmap
dir.create('figure/Pancan/')
pdf('figure/Pancan/Heatmap.pdf', width = 10, height = 10)
MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column = F,
                  anno_pch_cols = c('Histology', "Treatment"),
                  anno_regular = c("Institution","PIK3CA mutation","Passage (pro)", 'cluster'))
MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column_by = 'Histology')
MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column_by = 'cluster',
                  anno_pch_cols = c('Histology', "Treatment"),
                  anno_regular = c("Institution","PIK3CA mutation","Passage (pro)", 'cluster'))

MakeTopProHeatmap(pancan_obj, palette_list = palette_list_pancan, split_column_by = 'PIK3CA mutation',
                  anno_pch_cols = c('Histology', "Treatment"),
                  anno_regular = c("Institution","PIK3CA mutation","Passage (pro)", 'cluster'))
dev.off()
pdf('figure/Pancan/Umap.pdf', width = 10, height = 10)
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

