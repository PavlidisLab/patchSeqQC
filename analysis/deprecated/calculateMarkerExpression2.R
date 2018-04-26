source('R/calcCellContam.R')

### calculate marker sums in aibs data
print('Calculating marker expression for broad cell types in Tasic')

aibsExprDataDf$contam_type = aibsExprDataDf$norm_broad_type
aibs_contam_all_broad = calcContamAllTypes(aibsExprDataDf, fullMarkerList) 
aibs_contam_all_broad$contam_type = factor(aibsExprDataDf$contam_type)
aibs_med_exprs_broad = aibs_contam_all_broad %>% group_by(contam_type) %>% summarize_all(median) %>% as.data.frame()
rownames(aibs_med_exprs_broad) = aibs_med_exprs_broad$contam_type

aibsExprDataDf$contam_type = paste0(aibsExprDataDf$norm_sub_type, '_on')
aibs_contam_all_sub = calcContamAllTypes(aibsExprDataDf, fullMarkerList) 
aibs_contam_all_sub$contam_type = factor(aibsExprDataDf$contam_type)
aibs_med_exprs_sub = aibs_contam_all_sub %>% group_by(contam_type) %>% summarize_all(median) %>% as.data.frame()
rownames(aibs_med_exprs_sub) = aibs_med_exprs_sub$contam_type

aibs_med_exprs = rbind(aibs_med_exprs_broad, aibs_med_exprs_sub[c('Ndnf_on', 'Sncg_on', 'Pvalb_on', 'Pyramidal_on'), ])

# compare_cell_types = setdiff(names(fullMarkerList), c('Ndnf', 'Pvalb', 'Sncg')) %>% make.names()
# expected_expr = aibs_med_exprs['Ndnf', compare_cell_types] %>% unlist
# 
# compare_cell_type_exprs = aibs_med_exprs[compare_cell_types, compare_cell_types] %>% as.matrix %>% diag


#### calculate expected marker sums in Zeisel data
print('Calculating marker expression for broad cell types in Zeisel')


zeisel_expr = zeiselExprDataDf %>% as.data.frame()
colnames(zeisel_expr) = make.names(colnames(zeisel_expr))
zeisel_expr_norm = zeisel_expr
zeisel_expr_norm[, zeisel_gene_names %>% make.names] = zeisel_expr[, zeisel_gene_names %>% make.names] *1E6 / (zeisel_expr[, zeisel_gene_names %>% make.names] %>% rowSums )
zeisel_expr_norm[, zeisel_gene_names %>% make.names] = zeisel_expr_norm[, zeisel_gene_names %>% make.names] + 1

zeisel_expr_norm$contam_type = zeisel_expr_norm$norm_broad_type
zeisel_contam_all_broad = calcContamAllTypes(zeisel_expr_norm, fullMarkerList) 
zeisel_contam_all_broad$contam_type = factor(zeisel_expr_norm$contam_type)
zeisel_med_exprs_broad = zeisel_contam_all_broad %>% group_by(contam_type) %>% summarize_all(median) %>% as.data.frame()
rownames(zeisel_med_exprs_broad) = zeisel_med_exprs_broad$contam_type

zeisel_expr_norm$contam_type = paste0(zeisel_expr_norm$norm_sub_type, '_on')
zeisel_contam_all_sub = calcContamAllTypes(zeisel_expr_norm, fullMarkerList) 
zeisel_contam_all_sub$contam_type = factor(zeisel_expr_norm$contam_type)
zeisel_med_exprs_sub = zeisel_contam_all_sub %>% group_by(contam_type) %>% summarize_all(median) %>% as.data.frame()
rownames(zeisel_med_exprs_sub) = zeisel_med_exprs_sub$contam_type

zeisel_med_exprs = rbind(zeisel_med_exprs_broad, zeisel_med_exprs_sub[c('Ndnf_on', 'Sncg_on', 'Pvalb_on', 'Pyramidal_on'), ])



### save marker sum values to the joined df per patch seq dataset

### per patch seq dataset and each cell within, calculate sum of marker expression

print('Calculating marker expression for all cells in patch seq data')

protein_coding_genes = annot %>% filter(gene_biotype == 'protein_coding') %>% select(mgi_symbol) %>% unlist %>% make.names
mt_gene_names = annot %>% filter(chromosome_name == 'MT') %>% select(mgi_symbol) %>% unlist %>% make.names

compare_cell_types_inh = setdiff(names(fullMarkerList), c('Ndnf_on', 'Pvalb_on', 'Sncg_on', 'Pyramidal_on', 'Inhibitory')) %>% make.names()
compare_cell_types_exc = setdiff(names(fullMarkerList), c('Pyramidal_on', 'Pvalb_on', 'Sncg_on', 'Ndnf_on', 'Pyramidal')) %>% make.names()

exc_cell_types = c('CA1-PYR', 'Exc', 'RS-PYR', 'BS-PYR')

for (i in 1:length(patch_seq_datasets)){
  dataset = patch_seq_datasets[[i]]
  dataset$joined_df$marker_sum = 0
  
  dataset_cell_types = levels(dataset$joined_df$major_type)
  marker_sums = lapply(dataset_cell_types, function(cell_type){
    
    curr_marker_type = dataset$joined_df[dataset$joined_df$major_type == cell_type, 'contam_type'][1] %>% as.character
    curr_marker_type = paste0(curr_marker_type, '_on')
    curr_marker_list = fullMarkerList[[curr_marker_type]]
    
    print(curr_marker_type)

    cell_inds = which(dataset$joined_df$major_type == cell_type)
    
    df = dataset$joined_df[cell_inds, ]
    rownames(df) = df$cell_id
    
    marker_expr = sumExpression(df, curr_marker_list)
    
    if (cell_type %in% exc_cell_types){
      compare_cell_types = compare_cell_types_exc
    } else{
      compare_cell_types = compare_cell_types_inh
    }
    
    if (names(patch_seq_datasets)[i] == 'fuzik'){
      compare_cell_types = compare_cell_types[compare_cell_types != 'OPC']
      compare_expr_dataset = zeisel_med_exprs[names(zeisel_med_exprs) != 'OPC', names(zeisel_med_exprs) != 'OPC']
      compare_expr_profile = zeisel_expr_norm[zeisel_expr_norm$contam_type == curr_marker_type, 
                                            fullMarkerList[c(curr_marker_type, compare_cell_types)] %>% unlist %>% unique] %>% log2 %>% colMeans
    }else{
      compare_expr_dataset = aibs_med_exprs
      compare_expr_profile = aibsExprDataDf[aibsExprDataDf$contam_type == curr_marker_type, 
                                            fullMarkerList[c(curr_marker_type, compare_cell_types)] %>% unlist %>% unique] %>% log2 %>% colMeans
    }
    # print(compare_expr_profile)

    contam_values = calcContamAllTypes(df, fullMarkerList)
    contam_values$contam_type = curr_marker_type
    
    contam_sum = normalizeContamSum(contam_values, compare_expr_dataset, compare_cell_types)
    marker_sum_norm = normalizeContam(contam_values, compare_expr_dataset, c(curr_marker_type, compare_cell_types))
    # marker_sum_norm_vec = marker_sum_norm[curr_marker_type, ] + 1 # need to add 1 because marker isn't contamination
    marker_sum_norm_vec = contam_values[, curr_marker_type] / compare_expr_dataset[curr_marker_type, curr_marker_type]
    

    dissoc_corr = rbind(compare_expr_profile %>% t() %>% as.data.frame, df[, fullMarkerList[c(curr_marker_type, compare_cell_types)] %>% unlist %>% unique] %>% log2) %>% t() %>% cor(method = 'spearman')
    dissoc_corr = dissoc_corr[1, -1]

    qc_measures = calculateDatasetQCMeasures(df, dataset$gene_names)

    out_df = data.frame(marker_sum = marker_expr, marker_sum_norm = marker_sum_norm_vec, contam_sum = contam_sum, cell_id = rownames(df), 
                        dissoc_corr = dissoc_corr)
    out_df = cbind(out_df, marker_sum_norm %>% t())
    out_df = merge(out_df, qc_measures, by = "cell_id")
    
    return(out_df)
  })
  marker_sums = bind_rows(marker_sums)
  
  marker_sum_df = bind_cols(dataset$joined_df %>% select(one_of('cell_id')), marker_sums)
  patch_seq_datasets[[i]]$joined_df$marker_sum =marker_sum_df$marker_sum
  patch_seq_datasets[[i]]$joined_df$marker_sum_norm =marker_sum_df$marker_sum_norm
  patch_seq_datasets[[i]]$joined_df$contam_sum =marker_sum_df$contam_sum
  patch_seq_datasets[[i]]$joined_df$num_genes =marker_sum_df$num_genes
  patch_seq_datasets[[i]]$joined_df$mt_gene_expr =marker_sum_df$mt_gene_expr
  patch_seq_datasets[[i]]$joined_df$dissoc_corr =marker_sum_df$dissoc_corr
  
  patch_seq_datasets[[i]]$contam_scores = marker_sum_df# %>% select(one_of(c('cell_id', names(fullMarkerList))))
  rownames(patch_seq_datasets[[i]]$contam_scores) = dataset$joined_df$cell_id
  patch_seq_datasets[[i]]$contam_scores = marker_sum_df# %>% select(one_of(c('cell_id', names(fullMarkerList))))
  patch_seq_datasets[[i]]$contam_scores = merge(patch_seq_datasets[[i]]$contam_scores, dataset$joined_df[, c('cell_id', 'major_type')])
  # 
  # qc_measures = calculateDatasetQCMeasures(dataset$joined_df, dataset$gene_names)
  # 
  # patch_seq_datasets[[i]]$joined_df$num_genes =qc_measures$num_genes
  # patch_seq_datasets[[i]]$joined_df$mt_gene_expr =qc_measures$mt_gene_expr
}


