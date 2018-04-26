library(ppcor)
library(parallel)
library(wCorr)
library(weights)

prepareEphysMatrix = function(input_data_matrix, ephys_var_names, log_ephys_vars = log_vars, scale_columns = T){
  
  response_vars = ephys_var_names
  ephys.col.inds = which(colnames(input_data_matrix) %in% response_vars)
  ephys_mat = as.data.frame(input_data_matrix[,ephys.col.inds])
  
  for (k in 1:length(response_vars)){
    curr_y_var = response_vars[k]
    y_data = input_data_matrix[,curr_y_var]
    if (curr_y_var %in% log_vars){
      y_data = log10(y_data)
    }
    
    ephys_mat[,k] <- y_data
  }
  
  if (scale_columns == T){
    ephys_mat = scale(ephys_mat)
  }
  else{
    ephys_mat = ephys_mat
  }
  rownames(ephys_mat) = input_data_matrix$sample.name
  return(ephys_mat)
}

# write a function to prepare expression data matrix
prepareExprMatrix = function(input_data_matrix, gene_names, scale_columns = T, threshold_cols = T, expr_mean_threshold = 7.5, expr_sd_threshold = .75){
  
  x_data <- as.matrix(input_data_matrix[,gene_names])
  
  raw.expr.mat = x_data
  
  if (threshold_cols){
    expr.mat.sd = apply(raw.expr.mat, 2, sd, na.rm = T)
    expr.mat.mean = apply(raw.expr.mat, 2, mean, na.rm = T)      
    
    sd.thresh = expr_sd_threshold
    mean.thresh = expr_mean_threshold
    use.inds = expr.mat.mean > mean.thresh & expr.mat.sd > sd.thresh
    x_data = x_data[,use.inds]
  }
  
  if (scale_columns){
    expr_mat = scale(x_data)
  }
  else{
    expr_mat = x_data
  }
  rownames(expr_mat) = input_data_matrix$sample.name
  return(expr_mat)
}
# 
# calculateExprEphysCorr = function(expr_mat, ephys_mat, sample_count = 0, return_pvals, ephys_var_names){
#   response_vars = ephys_var_names
#   MIN_SAMPLE_COUNT = sample_count
#   
#   gene_corr_df_list = list()
#   for (i in 1:length(response_vars)){
#     curr_y_var = response_vars[i]
#     
#     ephys_vec = ephys_mat[,curr_y_var]
#     gene_corr_df = as.data.frame(matrix(nrow = ncol(expr_mat), ncol = 5))
#     colnames(gene_corr_df) = c("corr", "pval", "gene", "mn", "sd")
#     rownames(gene_corr_df) = t(colnames(expr_mat))
#     gene_corr_df['gene'] = rownames(gene_corr_df)
#     
#     gene_corr_df = mclapply(expr_mat %>% as.data.frame(),  function(vec){
#       expr_vec = vec
#       sample_count = sum(!is.na(ephys_vec + expr_vec))
#       if (sample_count < MIN_SAMPLE_COUNT - 1){
#         next
#       }
#       t = cor.test(ephys_vec, expr_vec, na.action = na.omit, method = "spearman")
#       out_vec = data.frame(pval = t$p.value, corr = t$estimate)
#       # out_vec['pval'] = t$p.value
#       # out_vec['corr']= t$estimate
#       return(out_vec)
#     }, mc.cores = 30)
#     gene_corr_df = do.call(rbind, gene_corr_df)
#     gene_corr_df$gene = rownames(gene_corr_df)
# #     for (k in 1:ncol(expr_mat)){
# #       expr_vec = expr_mat[,k]
# #       sample_count = sum(!is.na(ephys_vec + expr_vec))
# #       if (sample_count < MIN_SAMPLE_COUNT - 1){
# #         next
# #       }
# #       t = cor.test(ephys_vec, expr_vec, na.action = na.omit, method = "spearman")
# #       gene_corr_df[k,'pval']= t$p.value
# #       gene_corr_df[k,'corr']= t$estimate
# # #       gene_corr_df[k,'mn']= mean(as.matrix(j[gene_corr_df[k,'gene']]), na.rm = T)
# # #       gene_corr_df[k,'sd']= sd(as.matrix(j[gene_corr_df[k,'gene']]), na.rm = T)
# #     }
#     fdrvals = p.adjust(gene_corr_df$pval, method = "fdr")
#     gene_corr_df$fdr = -log10(fdrvals)
# 
#     gene_corr_df_list[[i]] = gene_corr_df
#   }
#   names(gene_corr_df_list) <- response_vars
#   
#   return(gene_corr_df_list)
# }
# 
# calculateExprEphysPartialCorr = function(expr_mat, ephys_mat, partial_mat, sample_count = 0, return_pvals, ephys_var_names){
#   response_vars = ephys_var_names
#   MIN_SAMPLE_COUNT = sample_count
#   
#   gene_corr_df_list = list()
#   for (i in 1:length(response_vars)){
#     curr_y_var = response_vars[i]
#     
#     ephys_vec = ephys_mat[,curr_y_var]
#     gene_corr_df = as.data.frame(matrix(nrow = ncol(expr_mat), ncol = 5))
#     colnames(gene_corr_df) = c("corr", "pval", "gene", "mn", "sd")
#     rownames(gene_corr_df) = t(colnames(expr_mat))
#     gene_corr_df['gene'] = rownames(gene_corr_df)
#     for (k in 1:ncol(expr_mat)){
#       expr_vec = expr_mat[,k]
#       sample_count = sum(!is.na(ephys_vec + expr_vec))
#       if (sample_count < MIN_SAMPLE_COUNT - 1){
#         next
#       }
#       use_inds = !is.na(ephys_vec) & !is.na(expr_vec)
#       t = pcor.test(ephys_vec[use_inds], expr_vec[use_inds], partial_mat[use_inds, ], method = 'spearman')
#       #t = cor.test(ephys_vec, expr_vec, na.action = na.omit, method = "spearman")
#       gene_corr_df[k,'pval']= t$p.value
#       gene_corr_df[k,'corr']= t$estimate
#       #       gene_corr_df[k,'mn']= mean(as.matrix(j[gene_corr_df[k,'gene']]), na.rm = T)
#       #       gene_corr_df[k,'sd']= sd(as.matrix(j[gene_corr_df[k,'gene']]), na.rm = T)
#     }
#     fdrvals = p.adjust(gene_corr_df$pval, method = "fdr")
#     gene_corr_df$fdr = -log10(fdrvals)
#     
#     gene_corr_df_list[[i]] = gene_corr_df
#   }
#   names(gene_corr_df_list) <- response_vars
#   
#   return(gene_corr_df_list)
# }

# if (ephys_name %in% log_vars){
#   y_formula_term = paste0("log10(", ephys_name, ")")
# } else{
#   y_formula_term = ephys_name
# }
# formula = paste(y_formula_term, ' ~ log10(', gene, ')')

calculateExprEphysWeightedCorr = function(expr_mat, ephys_mat, weights, ephys_var_names, log_ephys_vars){
  
  # only perform calculation for genes with greater than 0 sd
  expr_mat_sd = apply(expr_mat, 2, sd)
  use_expr_mat = expr_mat[, expr_mat_sd > 0]
  
  # log10 transform expression data
  use_expr_mat = log2(use_expr_mat)
  
  all_corrs = mclapply(ephys_var_names, function(ephys_var){
    ephys_col = ephys_mat[, ephys_var]
    
    if (ephys_var %in% log_ephys_vars){
      ephys_col = log10(ephys_col)
    }
    use_inds = !is.na(ephys_col)
    t = wtd.cor(use_expr_mat[use_inds, ], ephys_col[use_inds], weights[use_inds])
    t_return = t %>% as.data.frame %>% tibble::rownames_to_column(var = 'gene')
    t_return$fdr = p.adjust(t_return$p.value, method = 'BH')
    t_return$ephys_prop = ephys_var
    return(t_return)
  }, mc.cores = 40)
  all_corrs = bind_rows(all_corrs) %>% select(gene, ephys_prop, correlation, p.value, fdr) %>% rename(corr = correlation, pval = p.value)
  return(all_corrs)
}


# now aggregate lists of differentially expressed genes
organizeCorrList = function(gene_corr_df_list, fdr_thresh =.05, ephys_var_names){
  response_vars = ephys_var_names
  
  fdr_thresh_2 = 1
  
  sig_genes_df = as.data.frame(matrix(nrow = 0, ncol = 0))
  sig_genes_df_all = as.data.frame(matrix(nrow = 0, ncol = 0))
  
  for (i in 1:length(response_vars)){
    curr_y_var = response_vars[i]
    if (nrow(gene_corr_df_list[[i]]) == 0){
      next
    }
    frame_out = gene_corr_df_list[[i]] %>% 
      dplyr::select(gene, pval, corr, fdr) %>%
      arrange(desc(fdr)) %>%
      mutate(ephys_name = curr_y_var)
    sig_genes_df_all = rbind(sig_genes_df_all, frame_out)
  }
  
  sig_genes_df = sig_genes_df_all %>% filter(fdr > -log10(fdr_thresh)) %>%
    arrange(ephys_name, desc(fdr), desc(corr))
  
  num_unique_sig_genes = sig_genes_df %>% 
    #filter(nt_adj_pval < .1) %>% 
    dplyr::select(gene) %>% 
    unique %>% 
    nrow
  
  l = list("sig_genes_df_all" = sig_genes_df_all, "sig_genes_df" = sig_genes_df, "num_unique_sig_genes" = num_unique_sig_genes)
  return(l)
}

