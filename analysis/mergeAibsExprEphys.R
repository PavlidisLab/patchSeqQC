# merge aibs ephys and expression datasets at the level of cre-lines

# library(dae)
# source('R/ephysExprPlots.R')


## set up rules for splitting by layers or not

# indicates for which cre-lines to split by layers (based on if there's at least 10 cells with both upper and lower gene expression)
lines_split = c('Calb2', 'Cux2', 'Htr3a', 'Pvalb', 'Sst', 'Vip')


## add a field to aibsExprDataDf indicating layer splits or not

upper_regions_expr = c('L1', 'L1-L2/3', 'L2/3', 'upper')
lower_regions_expr = c('L4', 'L5', 'L6', 'L6a', 'L6b', 'lower')
aibsExprDataDf = aibsExprDataDf %>% mutate('layer_coarse' = case_when((.$dissection %in% upper_regions_expr) &
                                                     (.$sample.name %in% lines_split) ~ 'upper', 
                                                   (.$dissection %in% lower_regions_expr) &
                                                     (.$sample.name %in% lines_split) ~ 'lower',
                                                   TRUE ~ "all")) 

# first summarize expression vales (by cre line for now)

keep_cell_types = c('GABA-ergic Neuron', 'Glutamatergic Neuron')
#C130083M11Rik:mt_X57780
tasic_expr = aibsExprDataDf %>% filter(tdTomato == 'positive', pass_qc_checks == 'Y', broad_type %in% keep_cell_types) 
# tasic_expr = merge(tasic_expr, aibsColors, by = "sample.name")
tasic_expr = left_join(tasic_expr, expr_cell_colors, by = "primary_type")

# perform log2 based normalization prior to calculating means
tasic_expr[, aibs_gene_names] = tasic_expr[, aibs_gene_names] %>% log2

aibs_expr_filt_stats = tasic_expr %>%
  group_by(sample.name, layer_coarse) %>% select(Hcn3, broad_type) %>% 
  mutate(num_cells_expr = n(), 
         num_cells_exc = sum(broad_type == 'Glutamatergic Neuron')) %>% 
  mutate(prop_exc = num_cells_exc / num_cells_expr) %>%
  select(num_cells_expr, prop_exc) %>% distinct(.keep_all = T)


tasic_expr_cre = tasic_expr %>% group_by(sample.name, layer_coarse) %>% summarize_at(aibs_gene_names, mean)
# rownames(tasic_expr_cre) = tasic_expr_cre$sample.name

tasic_expr_cre[, aibs_gene_names] = 2^tasic_expr_cre[, aibs_gene_names]


### summarize ephys data by cre line

primary_vis_regions = c('VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b')
upper_regions = c('VISp1', 'VISp2/3')
lower_regions = c('VISp4', 'VISp5', 'VISp6a', 'VISp6b')

# adds a field to aibsEphys that indicates whether the brain region is upper or lower or all
aibsEphys = aibsEphys %>% filter(brain_region %in% primary_vis_regions) %>% 
  mutate('layer_coarse' = case_when((.$brain_region %in% upper_regions) & 
                                      (.$sample.name %in% lines_split) ~ 'upper' , 
                                      (.$brain_region %in% lower_regions) & 
                                      (.$sample.name %in% lines_split)  ~ 'lower',
                                     TRUE ~ "all"))
MIN_CELLS = 1
temp_ephys_names = response_vars[response_vars %in% colnames(aibsEphys)]

aibsEphysFilt = aibsEphys %>% group_by(sample.name, layer_coarse) %>%
  filter(reporter == 'positive')

aibs_ephys_cre = aibsEphysFilt %>% mutate(num_cells = n()) %>% filter(num_cells >= MIN_CELLS) %>% select_(.dots = temp_ephys_names) %>% summarize_all(funs(mean(., na.rm = T)))
# rownames(aibs_ephys_cre) = aibs_ephys_cre$sample.name

aibs_ephys_filt_stats = aibsEphys %>% 
  filter(reporter == 'positive') %>%
  group_by(sample.name, layer_coarse) %>% mutate(num_cells_ephys = n()) %>% filter(num_cells_ephys >= MIN_CELLS) %>% select(num_cells_ephys, colors) %>% distinct(.keep_all = T)


## join aibs ephys and expression data
aibs_cre_joined = merge(tasic_expr_cre, aibs_ephys_cre, by = c("sample.name", "layer_coarse"), all = F, unique = T)
# rownames(aibs_cre_joined) = aibs_cre_joined$sample.name
colnames(aibs_cre_joined) = make.names(colnames(aibs_cre_joined))
aibs_cre_joined = merge(aibs_cre_joined, aibsColors, by = "sample.name")

aibs_cre_joined = left_join(aibs_cre_joined, aibs_expr_filt_stats , by = c("sample.name", "layer_coarse")) 
aibs_cre_joined = left_join(aibs_cre_joined, aibs_ephys_filt_stats %>% select(-colors), by = c("sample.name", "layer_coarse")) 

sample_size_weights = apply(cbind(aibs_cre_joined$num_cells_ephys, aibs_cre_joined$num_cells_expr), 1, function(x) sqrt(harmonic.mean(x)))
aibs_cre_joined$sample_size_weights = sample_size_weights
aibs_cre_joined = aibs_cre_joined %>% mutate('broad_type' = case_when(.$prop_exc > .5 ~ 'Exc', TRUE ~ "Inh"))
aibs_cre_joined$broad_type_colors = 'red'
aibs_cre_joined$broad_type_colors[aibs_cre_joined$broad_type == 'Exc'] = 'turquoise'

# unite 
new_sample_names = aibs_cre_joined %>% select(one_of('sample.name', 'layer_coarse')) %>% unite(col = sample.name, sample.name, layer_coarse) %>% select(sample.name)
aibs_cre_joined$sample.name = factor(new_sample_names %>% unlist %>% as.character())

saveRDS(aibs_cre_joined, file = 'data/aibs_cre_expr_ephys.rda')



