library(tidyverse)
library(ComplexHeatmap)
library(patchwork)

##############################
## LOAD data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Data/Phospho/')
site_level = read_tsv('phosphosite_matrix-log2_ratios-MD_norm(1).tsv')
gene_level = read_tsv('phosphogene_matrix-log2_ratios-MD_norm(1).tsv')
peptide_level = read_tsv('phosphopeptide_matrix-log2_ratios-MD_norm(1).tsv')
meta = read_tsv('../Metadata/Proteomics_Meta_lodaing_merged_02222022.tsv')

##############################
## Process data
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/phospho/')

site_level$Phosphosite.Index %>% str_extract('_\\d_\\d$') %>% table
site_level$Phosphosite.Index %>% str_extract('_\\d_\\d_(?=[S,T,Y].+$)') %>% table


# 1. Break down index
# Phosphosite.Index: follows format below
# Database | UniprotID | GENE_SPECIES_#1_#2_#3_#4_Phosphosite(s)
# 
# 1.	first possible phospho site
# 2.	last possible phospho site
# 3.	number of possible phospho sites
# 4.	number of localized phospho sites (localization done by Luciphor2 successfully)

# Note will see some warnings since 1_0 and 2_0 doesn't have LocalizedPhosphosites
idx_col = 'Phosphosite.Index'
site_level_df = site_level %>% separate(col = idx_col, sep = '[|,_]', into = c('Database','UniprotID','Gene','Species',
                                                                                              'FirstPossiblePhosphoSite','SecondPossiblePhosphoSite',
                                                                                              'NumbPossiblePhosphoSite','NumbLocalizedPhosphoSite','LocalizedPhosphosites'), remove = F)
site_level_df[site_level_df=='None'] = NA # Replace None to NA
# 2. Keep 'Localized only"
site_level_Localized_df =site_level_df %>% filter(!is.na(LocalizedPhosphosites))
message(str_glue("Filtered out {nrow(site_level_df) - nrow(site_level_Localized_df)} fragments, with {nrow(site_level_Localized_df)} rows left."))

# 3. remove less relevant column
col_rm = c('Modifications','Flanking.Sequence.Phosphosites.6mer')
site_level_Localized_df = site_level_Localized_df %>% select(-all_of(col_rm))

# For Long table
id_cols = names(site_level_Localized_df) %>% str_subset('PA|pool', negate =T)
site_flattened_long = site_level_Localized_df %>% pivot_longer(cols = -id_cols,
                                                          names_to = 'Sample aliquot ID', 
                                                          values_to = 'Phospho_Site_level') %>% 
  mutate(Phospho_Site_level = as.numeric(Phospho_Site_level)) %>% 
  filter(!is.na(Phospho_Site_level)) # Remove Those without reads

# Add Long table + meta
site_long_meta = left_join(site_flattened_long, meta, by = 'Sample aliquot ID') 

## OUTPUT
dir.create('processed_table/06282022/')
write_tsv(site_level_Localized_df, 'processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed.tsv')
write_tsv(site_long_meta, 'processed_table/06282022/phosphosite_matrix-log2_ratios-MD_norm(1)_localized_processed_long_w_meta.tsv')

