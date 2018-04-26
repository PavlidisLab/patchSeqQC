library(weights)
library(reshape2)
library(tidyr)
library(parallel)
library(wCorr)

### save marker sum values to the joined df per patch seq dataset

patch_seq_datasets$cadwell$joined_df$contam_type = 'Ndnf'
patch_seq_datasets$cadwell$joined_df$contam_type = factor(patch_seq_datasets$cadwell$joined_df$contam_type)

patch_seq_datasets$foldy$joined_df[str_detect(patch_seq_datasets$foldy$joined_df$major_type, 'RS-INT'), 'contam_type'] = 'Sncg'
patch_seq_datasets$foldy$joined_df[str_detect(patch_seq_datasets$foldy$joined_df$major_type, 'FS-INT'), 'contam_type'] = 'Pvalb'
patch_seq_datasets$foldy$joined_df[str_detect(patch_seq_datasets$foldy$joined_df$major_type, 'CA1-PYR'), 'contam_type'] = 'Pyramidal'

patch_seq_datasets$foldy$joined_df$contam_type = factor(patch_seq_datasets$foldy$joined_df$contam_type)

patch_seq_datasets$fuzik$joined_df[str_detect(patch_seq_datasets$fuzik$joined_df$major_type, 'Exc'), 'contam_type'] = 'Pyramidal'
patch_seq_datasets$fuzik$joined_df[str_detect(patch_seq_datasets$fuzik$joined_df$major_type, 'Inh'), 'contam_type'] = 'Ndnf'

patch_seq_datasets$fuzik$joined_df$contam_type = factor(patch_seq_datasets$fuzik$joined_df$contam_type)

### per patch seq dataset and each cell within, calculate sum of marker expression
source('R/calcCellContam.R')

for (i in 1:length(patch_seq_datasets)){
  dataset = patch_seq_datasets[[i]]
  dataset$joined_df$marker_sum = 0
  
  dataset_cell_types = unique(dataset$joined_df$major_type)
  marker_sums = lapply(dataset_cell_types, function(cell_type){
    
    curr_marker_type = dataset$joined_df[dataset$joined_df$major_type == cell_type, 'contam_type'][1] %>% as.character
    curr_marker_list = fullMarkerList[[curr_marker_type]]
    
    cell_inds = which(dataset$joined_df$major_type == cell_type)
    
    df = dataset$joined_df[cell_inds, ]
    rownames(df) = df$cell_id
    
    marker_expr = sumExpression(df, curr_marker_list)
    out_df = data.frame(marker_sum = marker_expr, cell_id = rownames(df))
    return(out_df)
  })
  marker_sums = bind_rows(marker_sums)
  
  marker_sum_df = bind_cols(dataset$joined_df %>% select(one_of('cell_id')), marker_sums)
  patch_seq_datasets[[i]]$joined_df$marker_sum =marker_sum_df$marker_sum
}

### calculate correlations for patch-seq data with and without weighting by contamination

calculateWeightedEphysGeneCorrs = function(dataset, cell_type_name, ephys_name, test_gene_list){
  
  joined_df = dataset$joined_df[dataset$joined_df$major_type == cell_type_name & !is.na(dataset$joined_df[, ephys_name]), ]
  
  if (! (ephys_name %in% dataset$ephys_names)){
    return(w[F, ]) # update this
  }
  
  use_gene_list = intersect(dataset$gene_names, test_gene_list)
  # 
  
  weighted_results = mclapply(use_gene_list%>% make.names(), function(gene){
    weighted_corr = weightedCorr(log10(joined_df[, gene]), joined_df[, ephys_name], 
                 method = "Pearson", weights = (joined_df[, 'marker_sum']))
    return(weighted_corr)
  }, mc.cores = 40) 
  weighted_results = do.call(rbind, weighted_results)
  rownames(weighted_results) = use_gene_list
  colnames(weighted_results) = 'weighted_results'
  
  unweighted_results = mclapply(use_gene_list %>% make.names(), function(gene){
    unweighted_corr = cor(log10(joined_df[, gene]), joined_df[, ephys_name], 
                                 method = "pearson", use = 'p')
    return(unweighted_corr)
  } , mc.cores = 40)
  unweighted_results = do.call(rbind, unweighted_results)
  rownames(unweighted_results) = use_gene_list
  colnames(unweighted_results) = 'unweighted_results'
  
  w = cbind(weighted_results, unweighted_results) %>% as.data.frame
  w = tibble::rownames_to_column(w, var = "gene")
  return(w)
}

dataset_list = c('cadwell', 'fuzik', 'foldy')
analyze_ephys_props = c('aphw')

analyze_cell_types = c('eNGC', 'SBC', 'RS-INT', 'FS-INT', 'Inh')

weighted_result_df_all = lapply(patch_seq_datasets, function(dataset){
  measured_ephys_props = dataset$ephys_names
  valid_ephys_props = intersect(analyze_ephys_props, measured_ephys_props)
  use_cell_types = unique(dataset$joined_df$major_type)
  use_cell_types = intersect(use_cell_types, analyze_cell_types)
  
  weighted_results_list = lapply(valid_ephys_props, function(ephys_prop){
    r1 = lapply(use_cell_types, function(cell_type){
      calculateWeightedEphysGeneCorrs(dataset, cell_type, ephys_prop, test_gene_list)
    })
    names(r1) = use_cell_types
    weighted_result_df = bind_rows(r1, .id = "use_cell_types")
  })
  names(weighted_results_list) = valid_ephys_props
  weighted_result_df = bind_rows(weighted_results_list, .id = "ephys_prop")
})
names(weighted_result_df_all) = names(patch_seq_datasets)
weighted_result_df_all  = weighted_result_df_all %>% bind_rows(.id = "dataset") 

## calculate how highly expressed and variable each gene is in a patch-seq dataset

source('R/calculateGeneDatasetStats.R')

num_mean_bins = 10
num_var_bins = 10

use_gene_list = aibs_gene_names
stats_dfs = lapply(patch_seq_datasets, function(dataset){
  use_cell_types = unique(dataset$joined_df$major_type)
  use_cell_types = intersect(use_cell_types, analyze_cell_types)
  
  r1 = lapply(use_cell_types, function(cell_type){
    stats = calculateGeneDatasetStats(dataset, cell_type, intersect(dataset$gene_names, use_gene_list))
    
    stats[stats$mean_expr == 1, 'mean_expr'] = NA
    stats[stats$var_expr == 0, 'var_expr'] = NA
    
    stats = stats %>% filter(!is.na(mean_expr)) %>% mutate(mean_bin = ntile(mean_expr, num_mean_bins))
    s = lapply(1:num_mean_bins, function(bin_ind){
      return(stats %>% filter(mean_bin == bin_ind) %>% mutate(var_bin = ntile(var_expr, num_var_bins)))
      # return(stats %>% filter(mean_bin == bin_ind) %>% mutate(var_bin = percent_rank(mean_expr)))
    }) %>% bind_rows()
    
    return(s)
  })
  names(r1) = use_cell_types
  stat_df = bind_rows(r1, .id = "use_cell_types")
  
  return(stat_df)
})
names(stats_dfs) = names(patch_seq_datasets)
stats_dfs  = stats_dfs %>% bind_rows(.id = "dataset") 

# merge stats_dfs with results from patch-seq correlation analysis
weighted_result_df_all_merged = merge(weighted_result_df_all, stats_dfs)

aibs_ephys_corr_pre_merge = aibs_ephys_corrs %>% 
  filter(ephys_prop == 'aphw') %>% 
  select(gene, ephys_prop, corr) %>% 
  rename(aibs_corr = corr)

weighted_result_df_all_merged = merge(weighted_result_df_all_merged , aibs_ephys_corr_pre_merge) 
weighted_result_df_all_merged = weighted_result_df_all_merged %>% mutate(cons_weighted = aibs_corr * weighted_results > 0, cons_unweighted = aibs_corr * unweighted_results > 0)


# 
# merge( , aibs_ephys_corrs %>% filter(ephys_prop == 'aphw') %>% select(gene, ephys_prop, corr))
# 
# weighted_result_df_all = merge(merge(weighted_result_df_all, stats_dfs) , aibs_ephys_corrs %>% filter(ephys_prop == 'aphw') %>% select(gene, ephys_prop, corr))
# weighted_result_df_all = weighted_result_df_all %>% mutate(cons_weighted = aibs_inh * weighted_results > 0, cons_unweighted = aibs_inh * unweighted_results > 0)
# 
# 
# aibs_ephys_corrs$dataset = 'aibs'
# aibs_ephys_corrs$use_cell_types = 'inh'
# aibs_ephys_corrs$weighted_results = aibs_ephys_corrs$corr
# 
# r = bind_rows(weighted_result_df_all, aibs_ephys_corrs %>% filter(ephys_prop == 'aphw') %>% select(dataset, ephys_prop, use_cell_types, gene, weighted_results))

# calculate accuracy
dropout_rate_allowance = 1
var_bin_threshold = 4
mean_bin_threshold = 4

summary_accuracy_counts = weighted_result_df_all_merged %>% filter(gene %in% test_gene_list) %>% group_by(dataset, ephys_prop, use_cell_types) %>% 
  filter(dropout_rate <= dropout_rate_allowance, 
         #std_expr > mean_expr) %>% 
         var_bin > var_bin_threshold, mean_bin > mean_bin_threshold) %>% 
  summarize(weighted_accuracy = sum(cons_weighted, na.rm = T) / length(cons_weighted %>% na.omit), 
            unweighted_accuracy = sum(cons_unweighted, na.rm = T) / length(cons_unweighted %>% na.omit), 
            gene_count = length(cons_unweighted %>% na.omit)) 

summary_accuracy_counts_melted = summary_accuracy_counts %>% select(dataset, use_cell_types, ephys_prop, weighted_accuracy, unweighted_accuracy) %>% melt(value.name = 'accuracy')
summary_accuracy_counts_melted$dataset = factor(summary_accuracy_counts_melted$dataset, levels = c('foldy', 'fuzik', 'cadwell'))
summary_accuracy_counts_melted = summary_accuracy_counts_melted %>% unite(cell_type, dataset, use_cell_types, remove = F)

bar_cons_plot = ggplot(summary_accuracy_counts_melted %>% filter(ephys_prop %in% c('aphw', 'maxfreq')), aes(x = cell_type, y = accuracy)) + geom_hline(yintercept = .5, linetype = 2, color = 'grey') + 
  geom_bar(aes(fill = variable), position = 'dodge', stat = 'identity') +
  ylab('Overall consistency with AIBS') + facet_wrap(~ephys_prop)

bar_cons_plot


agg_result = weighted_result_df_all_merged %>% filter(ephys_prop == 'aphw') %>% filter(cons_weighted == T, mean_bin > mean_bin_threshold -2, var_bin > var_bin_threshold - 2, abs(weighted_results) > .01) %>% 
  group_by(gene) %>% summarize(count = n()) %>% filter(count > 2) %>% as.data.frame()


any_cons_gene_meeting_thresholds = weighted_result_df_all_merged %>% filter(ephys_prop == 'aphw') %>% 
  filter(cons_weighted == T, mean_bin > mean_bin_threshold, var_bin > var_bin_threshold, abs(weighted_results) > .01) %>% 
  distinct(gene) %>% unlist

aphw_dataset_genes = lapply(dataset_list, function(dataset_name){
  print(dataset_name)
  # weighted_result_df_all_merged %>% filter(ephys_prop == 'aphw', dataset == dataset_name) %>% filter(cons_weighted == T, mean_bin > mean_bin_threshold -2, var_bin > var_bin_threshold - 2, abs(weighted_results) > .01) %>%
  #   select(gene) %>% unlist
  weighted_result_df_all_merged %>% filter(ephys_prop == 'aphw', dataset == dataset_name) %>% filter(gene %in% any_cons_gene_meeting_thresholds, cons_weighted == T, abs(weighted_results) > .01) %>%
    select(gene) %>% unlist
})
names(aphw_dataset_genes) = dataset_list

intersecting_genes_3_datasets = intersect(aphw_dataset_genes[[3]], intersect(aphw_dataset_genes[[1]], aphw_dataset_genes[[2]]))



## plot just foldy RSint data for one example gene
levels(jdf$colors) = c(levels(jdf$colors), 'purple')
jdf[jdf$sample.name == 'Vip', 'colors'] = 'purple'

plotPatchSeqWeighting = function(gene, ephys_prop, aibs_data = aibs_cre_joined, datasets = patch_seq_datasets, 
                                 dataset_name = 'foldy', plot_cell_type = 'all'){
  

  ephys_label = getEphysNiceName(ephys_prop)

  p0 = ggplot(data = aibs_data[aibs_data$broad_type == 'Inh', ]
              , aes_string(x = gene, y = ephys_prop,  weight = 'sample_size_weights', color = 'colors', group = 'broad_type')) + 
    geom_smooth(method = 'lm', aes(group=1), se = FALSE, color="grey", alpha = .5, linetype = 1) + 
    geom_point(size = 4) + 
    scale_color_identity() + 
    scale_x_continuous(trans = 'log10') + annotation_logticks(sides = "b")  + 
    theme(legend.position="none") 
  p0 = p0 + xlab(paste(gene, '(TPM + 1)') ) + ylab(ephys_label) + ggtitle('AIBS inhibitory cell types,\n (pooled)')
  
  joined_df = datasets[[dataset_name]]$joined_df
  joined_df$weights = 1/(joined_df$contam_sum+.2)
  
  if (plot_cell_type == 'all'){
    use_cell_types = unique(joined_df$major_type)
  }
  
  p1 = ggplot(joined_df[joined_df$major_type %in% use_cell_types, ], 
              aes_string(x = gene, y = ephys_prop, color = 'colors')) +
    geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey')   +
    geom_point(size = 2) +
    scale_x_continuous(trans = 'log10') + annotation_logticks(sides = "b")  + 
    scale_color_identity() + 
    theme(legend.position="none")
  p1 = p1 + xlab(paste(gene, '(TPM + 1)') )  + ylab(ephys_label) + ggtitle(paste(dataset_name, use_cell_types, 'cells,\n (unweighted)'))
  
  p2 = ggplot(joined_df[joined_df$major_type %in% use_cell_types, ], 
              aes_string(x = gene, y = ephys_prop, color = 'colors', size = 'weights')) +
    geom_smooth(method = "lm", aes(weight = weights), se = F, linetype = 1, color = 'grey') +
    geom_point() + 
    scale_x_continuous(trans = 'log10') + annotation_logticks(sides = "b")  + 
    scale_color_identity() + 
    theme(legend.position="none")
  p2 = p2 + xlab(paste(gene, '(TPM + 1)') )  + ylab(ephys_label) + ggtitle(paste(dataset_name, use_cell_types, 'cells,\n (weighted)'))
  
  p = plot_grid(p0, p1, p2, nrow = 1)
  
  return(p)
}

