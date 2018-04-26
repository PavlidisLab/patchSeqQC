### export cell type mappings

zeisel_broad_primary_types = zeiselExprDataDf[, c('norm_broad_type', 'level2class')] %>% unique() %>% 
  filter(norm_broad_type %in% names(fullMarkerList)) %>%
  group_by(norm_broad_type) %>% 
  summarize(zeisel_primary_types = paste(level2class, collapse = ', '))
aibs_broad_primary_types = aibsExprDataDf[, c('norm_broad_types', 'primary_type')] %>% mutate(norm_broad_type = norm_broad_types) %>% 
  unique() %>%   filter(norm_broad_type %in% names(fullMarkerList)) %>%
  group_by(norm_broad_type) %>% 
  summarize(tasic_primary_types = paste(primary_type, collapse = ', '))

broad_type_mapping = left_join(aibs_broad_primary_types, zeisel_broad_primary_types)

write.csv(broad_type_mapping, 'analysis/tables/broad_cell_type_mapping.csv')