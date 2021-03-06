
#' Plot the summarized marker expression for dissociated cell and patch seq data
#'
#' @param patch_seq_contam_df contamination data frame for patch seq samples
#' @param dissoc_cell_contam_df contamination data frame for dissociated cell samples
#' @param on_type string for which cell type is the on type - all patch seq samples should be of this type
#' @param off_type string for which cell type is the off type
#'
#' @return plot showing marker expression for both datasets
#' @export
#'
#' @examples
plotPatchSeqDissocCellMarkerExpr = function(patch_seq_contam_df,
                                            dissoc_cell_contam_df,
                                            on_type,
                                            off_type){


x_min = 0
x_max = max(dissoc_cell_contam_df[, on_type], patch_seq_contam_df[, on_type]) + 0
y_min_1 = min(c(dissoc_cell_contam_df[, off_type], patch_seq_contam_df[,off_type], dissoc_cell_contam_df[,off_type])) - 0
y_max_1 = max(c(dissoc_cell_contam_df[,off_type], patch_seq_contam_df[,off_type], dissoc_cell_contam_df[,off_type])) + 0

PLOT_THRESHOLD = .025

aibs_pyr = dissoc_cell_contam_df %>% dplyr::filter(contam_type == off_type) %>% dplyr::select(off_type) %>% unlist
aibs_pyr_low = aibs_pyr %>% quantile(PLOT_THRESHOLD)
aibs_pyr_high = aibs_pyr %>% quantile(1-PLOT_THRESHOLD)
aibs_ndnf = dissoc_cell_contam_df %>% dplyr::filter(contam_type == on_type) %>% dplyr::select(on_type) %>% unlist
aibs_ndnf_low = aibs_ndnf %>% quantile(PLOT_THRESHOLD)
aibs_ndnf_high = aibs_ndnf %>% quantile(1-PLOT_THRESHOLD)

p1 = ggplot(dissoc_cell_contam_df %>% dplyr::filter(contam_type %in% c(on_type, off_type)), aes_string(x = on_type, y = off_type, color = 'colors')) +
  geom_hline(yintercept = aibs_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) +
  geom_vline(xintercept = aibs_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
  geom_point(alpha = 1) + scale_color_identity() +
  xlab(paste(on_type, 'markers (log2 TPM+1)')) + ylab(paste(off_type, 'markers (log2 TPM+1)')) +
  ggtitle('Reference cell data') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))


p2 = ggplot(patch_seq_contam_df, aes_string(x = on_type, y = off_type, color = 'colors')) +
  geom_hline(yintercept = aibs_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) +
  geom_vline(xintercept = aibs_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() +
  xlab(paste(on_type, 'markers (log2 TPM+1)')) + ylab(paste(off_type, 'markers (log2 TPM+1)')) +
  ggtitle('Patch seq cell data') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

cadwell_aibs_marker_comp = plot_grid(p1, p2, nrow = 1)
return(cadwell_aibs_marker_comp)

}
