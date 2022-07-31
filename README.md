# pik3ca_proteomics
pik3ca proteomics phospho processing and analysis script

- This report host data processing and analysis script for PIK3CA project.
- More info can be found in Wiki

## KS Heatmap plot
- Added 7/30/2022
- From `/proteo_phospho/analysis/EGFR_score_GBMpaper.R`
```R
##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/')
protein_flattened_wide = read_tsv('proteomic/processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed.tsv')
phos_flattened_wide = read_tsv('phospho/processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed.tsv')
meta = read_tsv('../Data/Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')

kinase_table = read_tsv('../Resource/Marker/kindase_substrate_GBM_paper_Fig3d.tsv')


##############################
## Source functions
walk(list.files('script/proteomics/src/', full.names = T), source)
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/proteo_phospho/')

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
pdf('figure/Kinase_Substrate_heatmap_EGFR.pdf', width =10 ,height = 12)
MakeKSHeatmap(protein_flattened_wide, phos_flattened_wide, meta_use, kinase_table = kinase_table, 
              kinase_of_interest = 'EGFR',
              na_cutoff_kinase = 1, na_cutoff_substrate = 0.5)
dev.off()
```
<img width="755" alt="image" src="https://user-images.githubusercontent.com/54045654/182013094-685bd831-1f52-4e00-93ed-7a51a0eca490.png">
