
getGoodMarkersFxn = function(on_cell_name, markers, dissoc_cell_dataset = 'tasic', MARKER_THRESHOLD = 10, MIN_CELLS_EXPRESSING_RATIO = .75, 
                          marker_direction = 'positive', use_subtype = T){
  # checks whether found markers are expressed at a minimum level
  # and whether markers are expressed in at least XX percent of cells
  
  if(dissoc_cell_dataset == 'tasic'){
    dissoc_cell_expr = aibsExprDataDf
    # MARKER_THRESHOLD = 10 + 1 # in TPM in Tasic
    # MIN_CELLS_EXPRESSING_RATIO = .75 # in ratio
  } else{
    dissoc_cell_expr = zeiselExprDataDf
    # MARKER_THRESHOLD = 1 # in UMI in Zeisel
    # MIN_CELLS_EXPRESSING_RATIO = .5 # in ratio
  }
  
  # MARKER_THRESHOLD = 10 + 1 # in TPM in Tasic
  # MIN_CELLS_EXPRESSING_RATIO = .75 # in ratio
  
  if(use_subtype == F){
    num_cells_of_type = sum(dissoc_cell_expr$norm_broad_type == on_cell_name)
    expr_mat = dissoc_cell_expr[dissoc_cell_expr$norm_broad_type == on_cell_name, markers]
  }else{
    num_cells_of_type = sum(dissoc_cell_expr$norm_sub_type == on_cell_name)
    expr_mat = dissoc_cell_expr[dissoc_cell_expr$norm_sub_type == on_cell_name, markers]
  }
  
  markers_matrix =  expr_mat %>% 
    summarize_all(funs(sum(. >= MARKER_THRESHOLD))) %>% t() %>% 
    as.data.frame() %>% tibble::rownames_to_column(var = 'gene') %>%
    rename(nnz_gene_count = V1)
  
  matching_markers = markers_matrix %>% filter(nnz_gene_count >= num_cells_of_type * MIN_CELLS_EXPRESSING_RATIO) %>% 
    select(gene) %>% unlist %>% as.character
  
  if (marker_direction == 'negative'){
    matching_markers = setdiff(markers, matching_markers)
  }
  return(matching_markers)
}

