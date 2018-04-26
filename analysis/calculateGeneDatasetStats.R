## calculate statistics like means and variances per gene per cell type per patch seq dataset
source('R/calculateGeneDatasetStats.R')

num_mean_bins = 10
num_var_bins = 10

test_gene_list = aibs_gene_names
stats_dfs = lapply(patch_seq_datasets, function(dataset){
  use_cell_types = unique(dataset$joined_df$major_type)
  use_cell_types = intersect(use_cell_types, analyze_cell_types)
  
  r1 = lapply(use_cell_types, function(cell_type){
    stats = calculateGeneDatasetStats(dataset, cell_type, intersect(dataset$gene_names, test_gene_list))
    
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


# calcNonZeroCnt = function(expr_vec, zero_expr = 0){
#   return (sum(expr_vec != zero_expr))
# }
