########
# Palette 
#########
passage_col = colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(10) %>% setNames(0:9)
passage_col = c(passage_col, 'TBD' = '#666666', 'Pooled' = '#E78AC3')
palette_meta = list(Institution = c(WU = "#FFFF99", JHU = "#FB9A99", `MDACC-Kopetz` = "#80B1D3", 
                                      MDACC = "#6A3D9A", Wistar = "#E7298A", Welm = "#B2DF8A"),
                      Treatment = c(Sapanisertib = "#1B9E77", Cabozantinib = "#A65628"),
                      Histology = c(RCC = "#E78AC3", CRC = "#B3DE69", `Lung Adenocarcinoma` = "#DECBE4", 
                                    melanoma = "#8DD3C7", `breast cancer` = "#FED9A6"),
                      "PIK3CA mutation" = c(WT = "#CAB2D6", D350G = "#F0027F", H1047R = "#386CB0", `N345K, H1047R` = "#66C2A5", 
                                            Pooled = "#D9D9D9", E542K = "#E41A1C", Q546R = "#FFD92F", E545K = "#A6761D", 
                                            N1044K = "#6A3D9A", E726K = "#E78AC3", `E542K, E726K` = "#B3DE69", 
                                            `E545K (Node met)` = "#CCEBC5"),
                      'Passage (pro)' = passage_col,
                      'cluster' = c(`1` = "#FFD92F", `2` = "#66C2A5", `3` = "#377EB8", `4` = "#FC8D62", 
                                    `5` = "#FCCDE5")
)