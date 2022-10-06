library(tidyverse)
#library(edgeR)
library(DEP)

########################################
# Load CPTAC proteomics 
########################################
cptac_path = '~/Library/CloudStorage/Box-Box/Data_freeze_1.0'

# Load Proteomics
pro_root = '~/Library/CloudStorage/Box-Box/Data_freeze_1.0/Proteomics/CDAP_formatted/'
pho_paths = list.files(pro_root, pattern = 'PHO', recursive = T)
pho_name = pho_paths %>% str_remove('\\/.+')
pho_list = map(pho_paths, function(path){
  read_tsv(str_glue('{pro_root}/{path}'))
}) %>% setNames(pho_name)

# Load Mutation
mut_root = '~/Library/CloudStorage/Box-Box/Data_freeze_1.0/Somatic_mutation/WXS/'
mut_paths = list.files(mut_root, recursive = T, 'maf$')
mut_pahts = c("BRCA/CPTAC2_BRCA_prospective.v1.4.somatic.variants.070918.maf", 
              "ccRCC/ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf", 
              "CO/CPTAC2_CO_prospective.v1.3.somatic.variants.031918.maf", 
              "GBM/archived/tindaisy_all_cases_filtered.v2.0.20190905.maf", 
              "HNSCC/HNSCC_WXS_111.withmutect.maf", 
              "LSCC/lscc-v3.0-somaticvariant-SW/WashU.SW.V3.LSCC.Somatic.042420.maf", 
              "LUAD/LUAD.Somatic.050919.mnp.annot.maf", "OV/CPTAC2_Prospective_OV.v1.5.hg38.somatic.2019-01-19.maf", 
              "PDA/77_PDA_bulk.withmutect.withcaller.maf", "UCEC/ucecc.somatic.consensus.gdc.wu.113018.maf"
)
mut_name  = mut_pahts %>% str_remove('\\/.+')
mut_list = map(mut_pahts, function(path){
  read_tsv(str_glue('{mut_root}/{path}'))
}) %>% setNames(mut_name)

##################
# Process data
##################
FilterNumericMtxByNARatio = function(table, na_threshold = 0.5, verbose =T){
  # Numeric part
  mtx = table %>% select(where(is.numeric))
  
  # Count NA
  row_keep = is.na(mtx) %>% apply(1, mean) %>% `<`(., na_threshold)
  
  # report
  if(verbose) message(str_glue('Threshold {na_threshold} Filtered out {sum(!row_keep)} rows and have {sum(row_keep)} rows left'))
  # Filter
  cbind(table[row_keep, ])
}
pho_list_filtered = map(pho_list, FilterNumericMtxByNARatio, na_threshold = 0.6)
pho_mtx_list = map(pho_list_filtered[1], ~mutate(., Phos = str_c(Gene, '_', Phosphosite)) %>% 
                     # avg 'same' phosphosite
                     group_by(Phos) %>% summarise(across(where(is.numeric), mean, na.rm=T)) %>%
                     remove_rownames() %>%
                     column_to_rownames('Phos') %>% select(where(is.numeric))
                     )
          
# Load Expression - Not needed??
# Load Mutation

###################
# Extract Mutation
###################
GetMutStatus = function(maf, mtx, gene_use = 'PIK3CA'){
  # Numeric names = sample ids
  mtx_sample_id = mtx %>% select(where(is.numeric)) %>% names
  
  # Maf
  mut_df = maf %>% filter(str_detect(Hugo_Symbol, gene_use)) %>% select(c('Hugo_Symbol','Tumor_Sample_Barcode','HGVSp_Short')) %>%
    group_by(Tumor_Sample_Barcode) %>% 
    mutate(
      Tumor_Sample_Barcode = str_remove(Tumor_Sample_Barcode, '_T'), 
      HGVSp_Short = str_remove(HGVSp_Short, 'p.')) %>% 
    summarize(HGVSp_Short = toString(HGVSp_Short)) %>% 
    mutate(MutType = ifelse(str_detect(HGVSp_Short, '543|542|1047'), 
                            'recurrent','rare')) 
  
  # result
  merge(data.frame(Sample = mtx_sample_id), 
        mut_df, by.x = 'Sample', by.y = 'Tumor_Sample_Barcode', all.x=T) %>% 
    mutate(MutType = ifelse(is.na(MutType), 'WT', MutType))
}

mut_meta_list = map2(mut_list[names(pho_list)], pho_list, ~GetMutStatus(.x, .y))

#Count Mute percentage
mut_meta_list %>% imap(~.x %>% mutate(CancerType = .y)) %>% 
  bind_rows %>% 
  count(CancerType, MutType) %>% 
  pivot_wider(id_cols = CancerType, names_from = MutType, values_from = n, values_fill = 0) %>% 
  rowwise() %>%
  mutate(Total = sum(c_across(where(is.numeric)), na.rm=T)) %>% 
  mutate(Mut_percent = (rare+recurrent) / Total) %>% 
  arrange(desc(Mut_percent))

##################
# Differentail Expression
#################
## DE-Pro and Phos
## DEP workflow

## BRCA
mtx_use   = pho_mtx_list$BRCA
mtx_use = mtx_use %>% rownames_to_column('ID') %>% mutate(name = str_remove(ID, '_S')) 
mut_group = mut_meta_list$BRCA %>% mutate(Mut_status = ifelse(MutType == 'WT','WT','Mut')) #%>% {setNames(.$Mut_status, .$Sample)} 
# Fill na/nan with 0
mtx_use[is.na(mtx_use)] = NA

# Generate a SummarizedExperiment object using an experimental design
experimental_design <- mut_group %>% select(Sample, Mut_status) %>% setNames(c('label', 'condition')) %>% mutate(replicate = 1:n())
columns_idx = grep('ID|name', colnames(mtx_use), invert=T)
data_se <- make_se(mtx_use, columns_idx, experimental_design) 

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)

# Normalize the data
#data_norm <- normalize_vsn(data_se)

# Visualize normalization by boxplots for all samples before and after normalization
#plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_se)

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_se, fun = "MinProb", q = 0.01)

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "WT")
# Test all possible comparisons of samples
#data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", show_row_names = FALSE)

# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "Mut_vs_WT", label_size = 2, add_names = TRUE)


##################
# Differentail Expression 
# Wilcoxon  
#################
FindMarkerWilcox = function(mtx, meta_data, sample_col = 'Sample', group.by = 'Mut_status'){
  # Extract identity
  ident_list = meta_data[,c(sample_col, group.by)] %>% split(.[[group.by]]) %>% map(~.[[sample_col]])
  # split table
  
}
tmp = FindMarkerWilcox(mtx_use, mut_group)

wilcox.test()


###
##%>% select(!!rlang::sym(sample_col), !!rlang::sym(group.by)) %>% 
#    setNames(c('sample', 'group'))
###