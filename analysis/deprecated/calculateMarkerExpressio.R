

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

compare_cell_types_inh = setdiff(names(fullMarkerList), c('Ndnf', 'Pvalb', 'Sncg')) %>% make.names()
compare_cell_types_exc = setdiff(names(fullMarkerList), c('Pyramidal', 'Pvalb', 'Sncg')) %>% make.names()

exc_cell_types = c('CA1-PYR', 'Exc')

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
    
    if (curr_marker_type %in% exc_cell_types){
      compare_cell_types = compare_cell_types_exc
    } else{
      compare_cell_types = compare_cell_types_inh
    }
    contam_values = calcContamAllTypes(df, fullMarkerList)
    contam_values$contam_type = curr_marker_type

    contam_sum = normalizeContamSum(contam_values, aibs_med_exprs, compare_cell_types)
    
    qc_measures = calculateDatasetQCMeasures(df, dataset$gene_names)
    
    out_df = data.frame(marker_sum = marker_expr, contam_sum = contam_sum, qc_measures, cell_id = rownames(df))
    
    
    return(out_df)
  })
  marker_sums = bind_rows(marker_sums)
  
  marker_sum_df = bind_cols(dataset$joined_df %>% select(one_of('cell_id')), marker_sums)
  patch_seq_datasets[[i]]$joined_df$marker_sum =marker_sum_df$marker_sum
  patch_seq_datasets[[i]]$joined_df$contam_sum =marker_sum_df$contam_sum
  patch_seq_datasets[[i]]$joined_df$num_genes =marker_sum_df$num_genes
  patch_seq_datasets[[i]]$joined_df$mt_gene_expr =marker_sum_df$mt_gene_expr
}

### calculate marker sums in aibs data

my_group = aibsExprDataDfNew$primary_type
my_group[str_detect(my_group, 'Micr')] = "Microglia"
my_group[str_detect(my_group, 'Astr')] = "Astrocyte"
my_group[str_detect(my_group, 'Oligo')] = "Oligodendrocyte"
my_group[str_detect(my_group, 'OPC')] = "OPC"
my_group[str_detect(my_group, 'Endo')] = "Endothelial"
my_group[str_detect(aibsExprDataDfNew$broad_type, 'Gluta')] = "Pyramidal"
my_group[str_detect(aibsExprDataDfNew$broad_type, 'GABA')] = "Inhibitory"
my_group[str_detect(aibsExprDataDfNew$primary_type, 'Sncg')] = "Sncg"
my_group[str_detect(aibsExprDataDfNew$primary_type, 'Ndnf')] = "Ndnf"
my_group[str_detect(aibsExprDataDfNew$primary_type, 'Pvalb')] = "Pvalb"



my_group = factor(my_group)

aibsExprDataDfNew$contam_type = my_group
aibs_major_types = levels(my_group)


aibs_contam_all = calcContamAllTypes(aibsExprDataDf, fullMarkerList) 
aibs_contam_all$contam_type = factor(my_group)

aibs_med_exprs = aibs_contam_all %>% group_by(contam_type) %>% summarize_all(median) %>% as.data.frame()
rownames(aibs_med_exprs) = aibs_med_exprs$contam_type

compare_cell_type_exprs = aibs_med_exprs[compare_cell_types, compare_cell_types] %>% as.matrix %>% diag

compare_cell_types = setdiff(names(fullMarkerList), c('Ndnf', 'Pvalb', 'Sncg')) %>% make.names()
expected_expr = aibs_med_exprs['Ndnf', compare_cell_types] %>% unlist

normalizeContam = function(contam_values, aibs_med_exprs, compare_cell_types){
  
  expected_expr = aibs_med_exprs[contam_values$contam_type[1], compare_cell_types] %>% unlist()
  print(contam_values)
  print(expected_expr)
  print(contam_values$contam_type[1])
  compare_cell_type_exprs = aibs_med_exprs[compare_cell_types, compare_cell_types] %>% as.matrix %>% diag
  
  normalizedContamValues = apply(contam_values[, compare_cell_types], 1, function(x) (x - expected_expr) / compare_cell_type_exprs)
  print(normalizedContamValues)
  return(normalizedContamValues)

}

normalizeContamSum = function(contam_values, aibs_med_exprs, compare_cell_types){
  normalizedContamValues = normalizeContam(contam_values, aibs_med_exprs, compare_cell_types) %>% repWZero() %>% colSums
  return(normalizedContamValues)
}

normalizeContamSum(foldy_contam_vip, aibs_med_exprs, compare_cell_types)

repWZero = function(df){
  ndf = df
  ndf[df < 0] = 0
  return(ndf)
}

calculateDatasetQCMeasures = function(joined_df, dataset_genes, protein_genes = protein_coding_genes, mito_genes = mt_gene_names, min_expr_val = 1){
  
  num_genes = rowSums(joined_df[, getValidGenes(protein_genes, dataset_genes)] > min_expr_val)

  mt_gene_expr = sumExpression(joined_df, mito_genes %>% make.names())

  out_df = data.frame(num_genes = num_genes, mt_gene_expr = mt_gene_expr)
  return(out_df)
}

protein_coding_genes = annot %>% filter(gene_biotype == 'protein_coding') %>% select(mgi_symbol) %>% unlist %>% make.names(., unique = T)
num_genes = rowSums(joined_df[, protein_coding_genes] > 1)

mt_gene_names = annot %>% filter(chromosome_name == 'MT') %>% select(mgi_symbol) %>% unlist
