
#' Title plot marker expression in single patch-seq cells vs average expression of markers in dissoc cells
#'
#' @param patch_seq_df combined gene expression and metadata data frame for patch seq samples
#' @param dissoc_cell_df combined gene expression and metadata data frame for dissociated cellsamples
#' @param contam_type specify the "on" cell type
#' @param plot_cell_sample_ids which sample_ids to plot
#' @param plot_marker_types which types of cells to plot marker expression for
#' @param markers list of lists with marker genes
#'
#' @return a ggplot object with marker expression for various cells
#' @export
#'
#' @examples
plotPatchSeqVsDissocMarkerExpr = function(patch_seq_df,
                                          dissoc_cell_df,
                                          contam_type,
                                          plot_cell_sample_ids,
                                          plot_marker_types,
                                          markers){

  ### illustrate calculation of quality scores using correlation with dissociated cell data along markers
  fullMarkerListNdnf = markers[plot_marker_types]
  marker_df_ndnf = do.call(rbind, lapply(seq_along(fullMarkerListNdnf), function(i){
    data.frame(CLUSTER=names(fullMarkerListNdnf)[i], fullMarkerListNdnf[[i]])
  }))
  colnames(marker_df_ndnf) = c('cell_type', 'gene')
  marker_df_ndnf %<>% dplyr::mutate_at(c('cell_type', 'gene'), dplyr::funs(as.character(.)))

  markers_collapsed = markers[plot_marker_types] %>% unlist %>% unique
  aibs_ndnf = dissoc_cell_df[, c('norm_sub_types', markers_collapsed)] %>% dplyr::filter(norm_sub_types == contam_type)

  # average expression profile of the typical Ndnf cell from AIBS
  aibs_ndnf_mean = aibs_ndnf %>% dplyr::select(-norm_sub_types) %>% log2 %>% colMeans
  aibs_ndnf_mean = 2^aibs_ndnf_mean

  cadwell_expr = cadwell_engc_df[, c('sample_id', markers_collapsed)] %>% dplyr::filter(sample_id %in% plot_cell_sample_ids)

  rownames(cadwell_expr) = cadwell_expr$sample_id
  cadwell_expr = cadwell_expr %>% dplyr::select(-sample_id)

  cadwell_aibs_df = cbind(t(cadwell_expr), aibs_ndnf_mean) %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::left_join(., marker_df_ndnf, by = 'gene') %>%
    dplyr::left_join(., cell_type_colors, by = 'cell_type')

  cadwell_aibs_df = cadwell_aibs_df %>% tidyr::gather(sample_id, expr_val, plot_cell_sample_ids)

  cadwell_aibs_df_plot= cadwell_aibs_df %>% dplyr::filter(sample_id %in% plot_cell_sample_ids)
  cadwell_aibs_df_plot$sample_id = factor(cadwell_aibs_df_plot$sample_id, levels = plot_cell_sample_ids)

  ndnf_marker_illustration_plot = cadwell_aibs_df_plot  %>% ggplot(aes(x = expr_val, aibs_ndnf_mean, color = colors)) +
    geom_smooth(method = "lm", se = F, linetype = 2, color = 'black') +
    geom_point() + scale_y_continuous(trans = "log2") + scale_x_continuous(trans = "log2") + geom_abline(slope = 1, intercept = 0) +
    scale_color_identity("Cell type markers", labels = cadwell_aibs_df_plot$cell_type, breaks = cadwell_aibs_df_plot$colors, guide = "legend") +
    xlab('Patch-seq expr (TPM + 1)') + ylab('Mean dissoc. cell expr (TPM + 1)') +
    facet_wrap(~sample_id, scales = 'free_y')
  return(ndnf_marker_illustration_plot)
}


cell_type_colors = c(Ndnf_on = "indianred1", Pyramidal = "turquoise", Pyramidal_on = "turquoise",
                     Microglia = "gray50", Sncg_on = "purple", Astrocyte = "green4", Pvalb = "red", Pvalb_on = "red",
                     Inhibitory = "red",
                     Endothelial = "brown",
                     Oligodendrocyte = "sandybrown",
                     OPC = 'darkorange4') %>% as.data.frame %>% tibble::rownames_to_column(var = "cell_type")
colnames(cell_type_colors) = c('cell_type', 'colors')
