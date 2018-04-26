library(parallel)
detach("package:dplyr", unload=TRUE)
library(dplyr)
source('R/corr_funcs.R')
source('R/ephysExprPlots.R') # try removing this dependency
source('R/calculateGeneDatasetStats.R')

# for all genes in cadwell, calculate corrs with all ephys props

MEAN_BIN_THRESHOLD = 3 # how high (in terms of bins, out of 10) does mean level of gene need to be expressed
MARKER_SUM_THRESHOLD = 0 # how high does cell's markers need to be expressed
MIN_CONTAMINATION_VALUE = .1 # if cell expresses contamination lower than this amount, set to this amount

merged_ephys_gene_corrs = lapply(names(patch_seq_datasets), function(dataset_name){
  
  curr_dataset = dataset_name
  joined_df = patch_seq_datasets[[curr_dataset]]$joined_df
  
  use_cell_ids = joined_df[, c('major_type', 'cell_id', 'num_genes', 'contam_sum', 'marker_sum_norm', 'rmp', 'aphw')] %>% 
    filter(marker_sum_norm > MARKER_SUM_THRESHOLD, !is.na(rmp) | !is.na(aphw)) %>% select(cell_id) %>% unlist
  
  use_dataset_inds = which(joined_df$cell_id %in% use_cell_ids)
  
  gene_names = patch_seq_datasets[[curr_dataset]]$gene_names
  ephys_names = patch_seq_datasets[[curr_dataset]]$ephys_names
  
  print(ephys_names)
  
  expr_mat = joined_df[use_dataset_inds, gene_names]
  ephys_mat = joined_df[use_dataset_inds, ephys_names]
  
  # weights = 1/(lapply(joined_df$contam_sum[use_dataset_inds], function(x) max(x, MIN_CONTAMINATION_VALUE)) %>% unlist)
  weights = lapply(joined_df$dissoc_corr[use_dataset_inds], function(x) max(x, MIN_CONTAMINATION_VALUE)) %>% unlist
  
  # weights = joined_df$marker_sum_norm[use_dataset_inds]
  
  # print(paste('Calculating mean and variance of expression values for', dataset_name, 'genes'))
  dataset_stats_dfs = calculateMatrixGeneStats(expr_mat, gene_names) 
  
  use_genes = dataset_stats_dfs %>% filter(mean_bin > MEAN_BIN_THRESHOLD) %>% select(gene) %>% unlist
  
  print(paste('Calculating ephys and expression correlations for', dataset_name, 'data - with weighting'))
  
  dataset_ephys_corrs_weighted = calculateExprEphysWeightedCorr(expr_mat[, use_genes], ephys_mat, weights = weights, 
                                                                ephys_var_names = ephys_names, 
                                                                log_ephys_vars = log_vars)
  

  print(paste('Calculating ephys and expression correlations for', dataset_name, 'data - without weighting'))
  
  expr_mat = joined_df[use_dataset_inds ,gene_names]
  ephys_mat = joined_df[use_dataset_inds , ephys_names]
  
  no_weights = rep(1, nrow(ephys_mat))
  
  dataset_ephys_corrs_unweighted = calculateExprEphysWeightedCorr(expr_mat[, use_genes], ephys_mat, weights = no_weights, 
                                                                  ephys_var_names = ephys_names, 
                                                                  log_ephys_vars = log_vars)
  
  ephys_corrs_comb = merge(dataset_ephys_corrs_weighted, dataset_ephys_corrs_unweighted, by = c('gene', 'ephys_prop'))
  ephys_corrs_comb = merge(ephys_corrs_comb, dataset_stats_dfs, by = c('gene'))
  ephys_corrs_comb$dataset = dataset_name
  
  # ephys_corrs_comb = bind_rows("weighted" = dataset_ephys_corrs_weighted, "unweighted" = dataset_ephys_corrs_unweighted, .id = "groups")
  return(ephys_corrs_comb)
}) %>% bind_rows()
