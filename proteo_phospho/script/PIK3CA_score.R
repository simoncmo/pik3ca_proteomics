library(tidyverse)
library(patchwork)

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
markers = read_tsv('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Resource/Marker/KEGG_PI3K-AKI_genelist.txt')
# markers_pik3ca = read_tsv('../Resource/Marker/PI3K_score.tsv')
# https://www.genome.jp/pathway/hsa04151
marker_cynthia = read_tsv('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Resource/Marker/CynthiaTNBC_markers.txt')

##############################
# KS table
kinase_table = read_tsv('../Resource/Marker/kindase_substrate_GBM_paper_Fig3d.tsv')
# Change to protein name
kinase_table = kinase_table %>% mutate(Kinase = str_replace(Kinase, 'PDGFRA', 'PGFRA'))
kinase_table = kinase_table %>% mutate(Kinase = str_replace(Kinase, 'PTPN11', 'PTN11'))


##############################
# Speicies filter
protein_flattened_wide = protein_flattened_wide %>% filter(Species == "HUMAN")
protein_long_meta = protein_long_meta %>% filter(Species == "HUMAN")
phos_flattened_wide = phos_flattened_wide %>% filter(Species == "HUMAN")
phos_long_meta = phos_long_meta %>% filter(Species == "HUMAN")

##############################
## Use PI3K and relevent gene:
# PK3CA, PK3CB, PK3CD, PK3C3, AKT1,2,3, MTOR, PTEN
# Make fake Kinase Substrate table to run script below first
# Set All effector as 'substrate'
# including all PIK3CA isoform. ONly PK3CA as kinase

kinase_table_PI3K = data.frame(Kinase = 'PK3CA',
                              Substrate = c('PK3CA','PK3CB','PK3CD','PK3C3','PK3CG','AKT1','AKT2',"AKT3",'MTOR','PTEN')
)

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

#source('../script/shared/function_make_ks_heatmap.R')
source('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/script/shared/function_make_prophos_heatmap.R')


##############################
# Plot
##############################

########### 
# MAKE COMBINED Pro Phos table
prophos_wide_df = MakeCombinedProPhosMtx(protein_flattened_wide, phos_flattened_wide, gene_use = c('PK3CA','PTEN'))
prophos_wide_df = FilterMtxAndFillNa(prophos_wide_df)
prophos_wide_df = SortMtxHclust(prophos_wide_df)

pdf('figure/PI3K_heatmap_v2.pdf', width =10 ,height = 8)
MakeProPhosHeatmap(protein_flattened_wide, phos_flattened_wide, 
  split_by_gene = T,
  gene_use = c('PK3CA','PK3CB','PK3CD','PK3C3','PK3CG','AKT1','AKT2',"AKT3",'MTOR','PTEN'),
meta_table = meta_use
)

MakeProPhosHeatmap(protein_flattened_wide, phos_flattened_wide, 
  split_by_gene = F,
  gene_use = c('PK3CA','PK3CB','PK3CD','PK3C3','PK3CG','AKT1','AKT2',"AKT3",'MTOR','PTEN'),
meta_table = meta_use
)
dev.off()



pdf('figure/TNBC_marker_heatmap.pdf', width =10 ,height = 12)
MakeProPhosHeatmap(protein_flattened_wide, phos_flattened_wide, 
  split_by_gene = T,
  gene_use = marker_cynthia$Uniprot_name,
meta_table = meta_use
)

MakeProPhosHeatmap(protein_flattened_wide, phos_flattened_wide, 
  split_by_gene = F,
  gene_use = marker_cynthia$Uniprot_name,
meta_table = meta_use
)
dev.off()

library(ComplexHeatmap)
set.seed(12345)
m = matrix(rnorm(100), 20)
ha = HeatmapAnnotation(foo = anno_boxplot(m, height = unit(4, "cm")), which = 'row')
Heatmap(m, right_annotation = ha)
