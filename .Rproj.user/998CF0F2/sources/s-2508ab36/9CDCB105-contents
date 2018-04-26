

calculateGeneDatasetStats = function(dataset, cell_type, test_gene_list){
  
  use_gene_list = intersect(dataset$gene_names, test_gene_list)
  expr_dataset = dataset$joined_df[dataset$joined_df$major_type == cell_type, use_gene_list]
  mean_expr = apply(expr_dataset, 2, function(x) mean(log2(x)))
  var_expr = apply(expr_dataset, 2, function(x) var(log2(x)))
  std_expr = apply(expr_dataset, 2, function(x) sd(log2(x)))
  dropout_rate = apply(expr_dataset, 2, calcDropoutRate, 1)
  
  stat_df = cbind(mean_expr, var_expr, std_expr, dropout_rate) %>% as.data.frame() %>% tibble::rownames_to_column(var = 'gene')
  return(stat_df)
}

calcDropoutRate = function(expr_vec, zero_expr = 0){
  return (sum(expr_vec == zero_expr) / length(expr_vec))
}


calculateMatrixGeneStats = function(joined_df, use_gene_list, num_mean_bins = 10, num_var_bins = 10){
  
  expr_dataset = joined_df[, use_gene_list] %>% as.matrix()
  mean_expr = apply(expr_dataset, 2, function(x) mean(log2(x)))
  var_expr = apply(expr_dataset, 2, function(x) var(log2(x)))
  std_expr = apply(expr_dataset, 2, function(x) sd(log2(x)))
  dropout_rate = apply(expr_dataset, 2, calcDropoutRate, 1)
  
  stat_df = data.frame(as.matrix(cbind(mean_expr, var_expr, std_expr, dropout_rate))) %>% tibble::rownames_to_column(var = 'gene')

  stat_df = stat_df %>% filter(mean_expr > 1) %>% mutate(mean_bin = ntile(mean_expr, num_mean_bins))
  stat_df = lapply(1:num_mean_bins, function(bin_ind){
    return(stat_df %>% filter(mean_bin == bin_ind) %>% mutate(var_bin = ntile(var_expr, num_var_bins)))
  }) %>% bind_rows()

  return(stat_df)
}

num_mean_bins = 10
num_var_bins = 10

test_gene_list = aibs_gene_names
stats_dfs = lapply(patch_seq_datasets, function(dataset){
  use_cell_types = unique(dataset$joined_df$major_type)
  # use_cell_types = intersect(use_cell_types, analyze_cell_types)
  
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