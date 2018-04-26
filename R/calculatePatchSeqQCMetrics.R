

#' Title calculate QC metrics in patch-seq data
#'
#' @param patch_seq_df combined patch seq dataset containing expression and metadata values. must include a column labelleled 'major_type' and 'contam_type'
#' @param comparison_dataset which comparison dataset to use? only accommodates aibs for now
#' @param markers list of lists with marker genes
#'
#' @return a data frame with summarized expression values per cell type and other qc fields
#' @export
#'
#' @examples
calculatePatchSeqQCMetrics = function(patch_seq_df, comparison_dataset = 'aibs', markers){
  # dataset$joined_df$marker_sum = 0 # fix me

  dataset_cell_types = unique(patch_seq_df$major_type) %>% as.character()
  marker_sums = lapply(dataset_cell_types, function(cell_type){

    curr_marker_type = patch_seq_df[patch_seq_df$major_type== cell_type, 'contam_type'][1] %>% as.character
    curr_marker_type = paste0(curr_marker_type, '_on')
    curr_marker_list = markers[[curr_marker_type]]

    cell_inds = which(patch_seq_df$major_type == cell_type)

    df = patch_seq_df[cell_inds, ]
    rownames(df) = df$sample_id

    marker_expr = sumExpression(df, curr_marker_list)

    if (cell_type %in% exc_cell_types){
      compare_cell_types = compare_cell_types_exc
    } else{
      compare_cell_types = compare_cell_types_inh
    }

    if (comparison_dataset == 'fuzik'){
      compare_cell_types = compare_cell_types[compare_cell_types != 'OPC']
      compare_expr_dataset = zeisel_med_exprs[names(zeisel_med_exprs) != 'OPC', names(zeisel_med_exprs) != 'OPC']
      compare_expr_profile = zeisel_expr_norm[zeisel_expr_norm$contam_type == curr_marker_type,
                                              markers[c(curr_marker_type, compare_cell_types)] %>% unlist %>% unique] %>% log2 %>% colMeans
    }else{
      compare_expr_dataset = aibs_med_exprs
      compare_expr_profile = aibsExprDataDf[aibsExprDataDf$contam_type == curr_marker_type,
                                            markers[c(curr_marker_type, compare_cell_types)] %>% unlist %>% unique] %>% log2 %>% colMeans
    }
    # print(compare_expr_profile)

    contam_values = calcContamAllTypes(df, markers)
    contam_values$contam_type = curr_marker_type

    contam_sum = normalizeContamSum(contam_values, compare_expr_dataset, compare_cell_types)
    marker_sum_norm = normalizeContam(contam_values, compare_expr_dataset, c(curr_marker_type, compare_cell_types))
    # marker_sum_norm_vec = marker_sum_norm[curr_marker_type, ] + 1 # need to add 1 because marker isn't contamination
    marker_sum_norm_vec = contam_values[, curr_marker_type] / compare_expr_dataset[curr_marker_type, curr_marker_type]


    quality_score = rbind(compare_expr_profile %>% t() %>% as.data.frame, df[, markers[c(curr_marker_type, compare_cell_types)] %>% unlist %>% unique] %>% log2) %>% t() %>% cor(method = 'spearman')
    quality_score = quality_score[1, -1]

    # qc_measures = calculateDatasetQCMeasures(df, dataset$gene_names)

    out_df = data.frame(marker_sum = marker_expr, marker_sum_norm = marker_sum_norm_vec, contam_sum = contam_sum, sample_id = rownames(df),
                        quality_score = quality_score)
    out_df = cbind(out_df, marker_sum_norm %>% t())
    # out_df = merge(out_df, qc_measures, by = "sample_id")

    return(out_df)
  })
  marker_sums = dplyr::bind_rows(marker_sums) %>% dplyr::select(sample_id, dplyr::everything())
  marker_sums = merge(patch_seq_df %>% dplyr::select(dplyr::one_of('sample_id', 'major_type', 'contam_type')), marker_sums, by= 'sample_id')

  marker_sum_df = marker_sums
  marker_sum_df = marker_sum_df[order(match(marker_sum_df$sample_id,patch_seq_df$sample_id)), ]
  return(marker_sum_df)

}
