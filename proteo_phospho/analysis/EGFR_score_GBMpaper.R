library(tidyverse)
library(patchwork)
#library(DEqMS)
#library(fgsea)

library(circlize)
library(ComplexHeatmap)

##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/')
protein_flattened_wide = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed.tsv')
protein_long_meta = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed_long_w_meta.tsv')
phos_flattened_wide = read_tsv('phospho/processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed.tsv')
phos_long_meta = read_tsv('phospho/processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed_long_w_meta.tsv')
meta = read_tsv('../Data/Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')
# markers = read_tsv('../Resource/Marker/KEGG_PI3K-AKI_genelist.txt')
# markers_pik3ca = read_tsv('../Resource/Marker/PI3K_score.tsv')
# https://www.genome.jp/pathway/hsa04151

kinase_table = read_tsv('../Resource/Marker/kindase_substrate_GBM_paper_Fig3d.tsv')


##############################
## Source functions
walk(list.files('script/proteomics/src/', full.names = T), source)
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/proteo_phospho/')


##############################
# Useful parameters
sample_ids = protein_flattened_wide %>% names %>%  str_subset('PA')
meta_cols_selected = c('Histology', 'Institution','Passage (pro)','Treatment', 'PIK3CA mutation')
########
# Palette 
#########
col_anno = c("Institution", "Passage (pro)", "Treatment", 'Histology', "PIK3CA mutation")
palette_list_use = MakeMetaPalette(meta = meta, columns_use = col_anno, existing_palette_list = palette_meta)
palette_list_pancan = c(palette_list_use)

##############################
# Process and simplify Meta
meta_use = meta %>% 
  filter(`Sample aliquot ID`!= 'reference pool') %>% 
  column_to_rownames('Sample aliquot ID') %>% .[sample_ids,meta_cols_selected]

source('../script/shared/function_make_ks_heatmap.R')


##############################
# Plot
##############################
# Change to protein name
kinase_table = kinase_table %>% mutate(Kinase = str_replace(Kinase, 'PDGFRA', 'PGFRA'))
kinase_table = kinase_table %>% mutate(Kinase = str_replace(Kinase, 'PTPN11', 'PTN11'))


pdf('figure/Kinase_Substrate_heatmap_EGFR.pdf', width =10 ,height = 12)
MakeKSHeatmap(protein_flattened_wide, phos_flattened_wide, meta_use, kinase_table = kinase_table, 
              kinase_of_interest = 'EGFR',
              na_cutoff_kinase = 1, na_cutoff_substrate = 0.5)
dev.off()

pdf('figure/Kinase_Substrate_heatmap_PTN11.pdf', width =10 ,height = 10)
MakeKSHeatmap(protein_flattened_wide, phos_flattened_wide, meta_use, kinase_table = kinase_table, 
              kinase_of_interest = 'PTN11',
              na_cutoff_kinase = 1, na_cutoff_substrate = 0.5)
dev.off()











