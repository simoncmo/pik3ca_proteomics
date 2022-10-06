library(tidyverse)
library(ComplexHeatmap)
library(patchwork)

##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/proteomic/')
protein_flattened_wide = read_tsv('processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed.tsv')
protein_long_meta = read_tsv('processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed_long_w_meta.tsv')
meta = read_tsv('../../Data/Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')

##########################################
# PLOT
# Bot plot
p_bar_histology = protein_long_meta %>% filter(!is.na(Histology)) %>% 
  ggplot(aes(x = Histology, y = Protein_level, fill = Read_species)) + 
  geom_boxplot( color = 'gray40') + 
  cowplot::theme_cowplot() + 
  Seurat::RotatedAxis() 

# Bot plot - By Species
p_bar_speicies = protein_long_meta %>% filter(!is.na(Histology)) %>% 
  ggplot(aes(x = Read_species, y = Protein_level, fill = Histology)) + 
  geom_boxplot( color = 'gray40') + 
  cowplot::theme_cowplot() + 
  Seurat::RotatedAxis() 

# Output
pdf('figure/QC_protein_level_barplot.pdf', width = 15, height= 8)
(p_bar_histology | p_bar_speicies) + plot_annotation(title =  'QC - Protein level Speicies vs Histology')
dev.off()

########################################
## PICK3CA activation score
#######################################
pathway_table = read_tsv('../../Resource/pathway/PI3K_pathway_protein_phoso.tsv')

pathway_table %>% mutate(Exist = ifelse(Protein %in% protein_long_meta$Protein, 'Exist','No')) %>% 
  count(Pathway, Exist) %>% 
  pivot_wider(id_cols = Pathway, names_from = 'Exist', values_from = 'n')


#### Consensus?
      
        