library(tidyverse)
library(Seurat)
library(ggthemes)

today_date <- Sys.Date() %>% format("%m%d%Y")
setwd("~/Box/Ding_Lab/Projects_Current/2021/PIK3CA_PDX/Analysis/")
# Load table
pdxmeta <- read_tsv("/Users/simonmo/Box/PDX-Pilot/DataFreeze/v5.data_freeze/1.sampleInfo/sampleInfo.washU_b1-b8.pdmr.other.passed.v5.extra.20210102.tsv")
meta <- read_tsv("../data/table/PIK3CA_cohort_meta_11012021.tsv")

# Get PIK3CA Model id and select samples from PDX Pilot table
modelids <- paste(meta$ModelID, collapse = "|")
metapik <- pdxmeta %>% filter(
  Original_ModelID %in% meta$ModelID,
  !str_detect(Group, "Normal")
)

# Add this 743489_274_T
meta_tmp <- pdxmeta %>% filter(ModelID_proper_v3 == "743489_274_T")
meta_tmp <- meta_tmp %>% mutate(Original_ModelID = str_replace_all(ModelID_proper_v3, "_", "-"))
metapik <- rbind(metapik, meta_tmp)

metapik <- merge(data.frame(
  Original_ModelID = meta$ModelID,
  Histology = meta$Histology,
  PIK3CA_SampleID = meta[["New SampleID"]],
  PIK3CA_ModelID = meta[["New ModelID"]]
),
by = "Original_ModelID",
metapik, all = T
)

write_tsv(metapik, "../data/table/PIK3CA_cohort_PDXPilotMeta_11022021.tsv")

# Fill NA
metapik[is.na(metapik)] <- "NoData"
col_select <- c("Original_ModelID", "PIK3CA_SampleID", "PIK3CA_ModelID", "CancerType", "Gender", "Histology")
metapik_clean <- metapik %>%
  select(all_of(col_select)) %>%
  distinct()
metact <- metapik %>%
  count(Original_ModelID, DataType, Group) %>%
  pivot_wider(values_from = n, names_from = Group)
metact <- merge(metact, metapik_clean, by = "Original_ModelID", all = T)
# mutate(PDX = ifelse(is.na(PDX),0,PDX))

# Plot
p <- metact %>% 
  ggplot(aes(x = Original_ModelID, y = DataType)) +
  facet_grid(. ~ Histology, scales = "free", space = "free") +
  geom_point(aes(fill = !is.na(Human_Tumor),
                 shape = is.na(PDX)), shape = 21, color = "#59c285", size = 6, stroke = 1.5) +
  # geom_point(aes(size = Human_Tumor*3), shape=22,  fill=NA, color = '#4287f5', stroke = 2)+
  
  scale_fill_manual(values = c("TRUE" = "#f58d4c", "FALSE" = "white")) +
  geom_text(aes(label = PDX)) +
  theme_clean() +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0, angle = 45))
p
ggsave(plot = p, str_glue("../figures/Data_availability_{today_date}.pdf"), width = 13, height = 2.2)






### Modified in google sheet
### LOAD from full availablilty table
avail_full = read_tsv('../data/table/PIK3CA_cohort_googleSheet_11102021.tsv')

# + meta and output
avail_volumn = full_join(avail_full , 
          meta[c(3,4,5,6,28,29)] %>% dplyr::rename(Original_ModelID = ModelID)
          ) %>% 
  select(Histology, Institution, everything())

# output
write_tsv(avail_volumn, '../data/table/PIK3CA_cohort_googleSheet_volumn_11102021.tsv', na='')

## Make PDX only table to make it simpler
avail_PDXonly = avail_volumn %>% 
  select(-Human_Tumor) %>%
  pivot_wider(names_from = DataType, values_from = PDX) %>% 
  mutate(Missing = case_when(is.na(`RNA-Seq`)&is.na(WES)~'Both',
                             is.na(`RNA-Seq`)&!is.na(WES)~'RNA-Seq',
                             !is.na(`RNA-Seq`)&is.na(WES)~'WES',
                             T~''
                             )) %>% 
  select(c(names(avail_volumn)[1:3], 'WES','RNA-Seq','Missing'), everything())

write_tsv(avail_PDXonly, '../data/table/PIK3CA_cohort_googleSheet_pdxonly_11102021.tsv', na='')

