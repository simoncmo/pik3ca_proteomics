library(tidyverse)
library(ComplexHeatmap)
library(patchwork)

##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Data/Proteomics/')
gene_level = read_tsv('gene_matrix-log2_ratios-MD_norm(2).tsv')
protein_level = read_tsv('protein_matrix-log2_ratios-MD_norm(2).tsv')
peptide_level = read_tsv('peptide_matrix-log2_ratios-MD_norm(2).tsv')
meta = read_tsv('../Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')

##############################
## Process data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/proteomic/')

# 1. Count speicies 
protein_level_df = protein_level %>% mutate(Read_species = case_when(str_detect(Protein.Group.Accessions, 'MOUSE') & str_detect(Protein.Group.Accessions, 'HUMAN')~ 'HUMAN+MOUSE',
                                                                str_detect(Protein.Group.Accessions, 'MOUSE') ~'MOUSE',
                                                                str_detect(Protein.Group.Accessions, 'HUMAN')~ 'HUMAN',
                                                                str_detect(Protein.Group.Accessions, 'PIG')~ 'PIG'))

## Protein level - Note each row have multiple possible gene target 
# ; mean there are Mouse and Human both - Speices Column will keep record of this
protein_flattened_wide = protein_level_df %>% 
    separate_rows(Protein.Group.Accessions, sep = ';') %>%  # Split row by ;
    separate(col = 'Protein.Group.Accessions', into = c('Database','Uniprot_id','Protein_Speicies'), sep = '\\|', remove = F) %>% # Split column by | 
    separate(col = 'Protein_Speicies', into = c('Protein','Species'), sep = '_', remove = T)

# Found 1 PIG gene - filter out
protein_flattened_wide = protein_flattened_wide %>% filter(Species != 'PIG')

# Turn 'None' to NA
protein_flattened_wide[protein_flattened_wide=='None'] = NA

# For Long table
id_cols = c('Protein.Group.Accessions', 'Protein.Group.Names', 'Species', 'Database', 'Uniprot_id', 'Protein', 'Read_species')
protein_flattened_long = protein_flattened_wide %>% pivot_longer(cols = -id_cols,
                                                          names_to = 'Sample aliquot ID', 
                                                          values_to = 'Protein_level') %>% 
  mutate(Protein_level = as.numeric(Protein_level)) %>% 
  filter(!is.na(Protein_level)) # Remove Those without reads

# Add Long table + meta
protein_long_meta = left_join(protein_flattened_long, meta, by = 'Sample aliquot ID') 

## OUTPUT
dir.create('processed_table/06072022/')
write_tsv(protein_flattened_wide, 'processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed.tsv')
write_tsv(protein_long_meta, 'processed_table/06072022/protein_matrix-log2_ratios-MD_norm(2)_processed_long_w_meta.tsv')

