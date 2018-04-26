library(matrixStats)
# calculate contamination by looking at each types markers

sumMarkerExpression = function(expr_matrix, target_markers_name, comparison_markers_name){
  
  # calc sum expr for target markers
  target_marker_genes = getValidGenes(mouseMarkerGenes$Cortex[[target_markers_name]], colnames(expr_matrix))
  contam_marker_genes = getValidGenes(mouseMarkerGenes$Cortex[[comparison_markers_name]], colnames(expr_matrix))
  
  target_expr = rowSums(expr_matrix[, target_marker_genes])
  comp_expr = rowSums(expr_matrix[, contam_marker_genes])
  
  return(cbind(target_expr, comp_expr))
}

sumExpression = function(expr_matrix, markers, colwise = F){
  # returns a vector of summed expression values per cell

  if (!(colwise)){
    target_marker_genes = getValidGenes(markers, colnames(expr_matrix))
    target_expr = rowSums(expr_matrix[, target_marker_genes %>% make.names()] %>% log2)
  } else{
    target_marker_genes = getValidGenes(markers, rownames(expr_matrix))
    target_expr = colSums(expr_matrix[target_marker_genes, ] %>% log2)
  }
  
  return(target_expr)
  
}

getValidGenes = function(input_genes, target_genes = aibs_gene_names){
  # returns set intersection of genes
  
  matching_inds = which(input_genes %in% target_genes)
  
  return(input_genes[matching_inds])
}

calcContamAllTypes = function(expr_matrix, markerList, colwise = F){
  # calculate sum of expression values given a list of marker genes per cell type
  
  contam_vals = lapply(markerList, function(cell_type_markers){
    return(sumExpression(expr_matrix, cell_type_markers, colwise))
  }) %>% as.data.frame()
  return(contam_vals)

}


normalizeContam = function(contam_values, aibs_med_exprs, compare_cell_types){
  # normalize cell's contamination values against expected values from tasic or comparison dataset
  
  expected_expr = aibs_med_exprs[contam_values$contam_type[1], compare_cell_types] %>% unlist()
  
  compare_cell_type_exprs = aibs_med_exprs[compare_cell_types, compare_cell_types] %>% as.matrix %>% diag
  
  normalizedContamValues = apply(contam_values[, compare_cell_types], 1, function(x) (x - expected_expr) / (compare_cell_type_exprs- expected_expr))
  return(normalizedContamValues)
  
}

normalizeContamSum = function(contam_values, aibs_med_exprs, compare_cell_types){
  # calculate sum of contamination values across all cell types in compare_cell_types
  normalizedContamValues = normalizeContam(contam_values, aibs_med_exprs, compare_cell_types) %>% repWZero() %>% colSums
  return(normalizedContamValues)
}


repWZero = function(df){
  # replace negative values with 0
  
  ndf = df
  ndf[df < 0] = 0
  return(ndf)
}

calculateDatasetQCMeasures = function(joined_df, dataset_genes, protein_genes = protein_coding_genes, mito_genes = mt_gene_names, min_expr_val = 1){
  # given an input patch seq data frame, and set of genes, calculate quality control measures
  
  valid_genes = getValidGenes(protein_genes, dataset_genes)
  num_genes = rowSums(joined_df[, valid_genes ] > min_expr_val)

  mt_gene_expr = sumExpression(joined_df, mito_genes %>% make.names())
  
  out_df = data.frame(num_genes = num_genes, mt_gene_expr = mt_gene_expr, cell_id = joined_df$cell_id)
  return(out_df)
}

