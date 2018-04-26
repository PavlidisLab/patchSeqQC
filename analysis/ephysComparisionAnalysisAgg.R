source('R/ephysExprPlots.R') # try removing this dependency

### Calculate how many genes are significantly associated with each ephys property in each dataset

FDR_THRESHOLD = .1
MEAN_BIN_THRESHOLD = 3

patchseq_sig_gene_counts = merged_ephys_gene_corrs %>% group_by(dataset, ephys_prop) %>% 
  summarize(weighted = sum(fdr.x < FDR_THRESHOLD & mean_bin > MEAN_BIN_THRESHOLD), 
            unweighted = sum(fdr.y < FDR_THRESHOLD & mean_bin > MEAN_BIN_THRESHOLD), 
            total_sig_count = weighted + unweighted) %>% filter(total_sig_count > 0)

patchseq_sig_gene_counts = patchseq_sig_gene_counts %>% select(-total_sig_count)%>% gather(key = c(dataset, ephys_prop), count, weighted:unweighted)
colnames(patchseq_sig_gene_counts)[3] = 'weighted_val'


response_vars_plot_names = c('Rin', 'Vrest', 'APthr', 'APamp', 'APhw', 'Tau', 'AHPamp', 'Rheo', 'FRmax', 'Cm', 'SFA')
patchseq_sig_gene_counts$ephys_prop = plyr::mapvalues(patchseq_sig_gene_counts$ephys_prop, from = response_vars, to = response_vars_plot_names)

dataset_formal_names = c('Cadwell', 'Földy', 'Fuzik', 'Tasic/AIBS')
dataset_unformal_names = c('cadwell', 'foldy', 'fuzik', 'aibs')
patchseq_sig_gene_counts$dataset = plyr::mapvalues(patchseq_sig_gene_counts$dataset, from = dataset_unformal_names, to = dataset_formal_names)

# set max bar
# patchseq_sig_gene_counts$count[patchseq_sig_gene_counts$count > 101] = 101

# plot histogram showing number of significant genes per ephys prop per dataset
patch_seq_gene_count_fig = patchseq_sig_gene_counts %>% ggplot(aes(x = factor(ephys_prop), y = count, fill = factor(weighted_val))) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  ylab('Gene count (FDR < 0.1)') + xlab('Ephys Property') + 
  scale_fill_manual("", values = c("unweighted" = "grey50", "weighted" = "black")) +
  facet_grid(~dataset, scales = 'free_x', space = 'free_x') + theme(legend.position = "none")
patch_seq_gene_count_fig


# go through significant genes per dataset and ephys property, count up overlaps and calculate p-values
# for expected overlaps

merged_ephys_gene_corrs_w_aibs = merge(merged_ephys_gene_corrs, aibs_ephys_corrs, by = c('gene', 'ephys_prop'))

merged_ephys_gene_corrs_w_aibs = merged_ephys_gene_corrs_w_aibs %>% mutate(cons_weighted = corr.x * corr > 0,
                                                                           cons_unweighted = corr.y * corr > 0)

patch_seq_tasic_consistency = merged_ephys_gene_corrs_w_aibs %>% group_by(dataset, ephys_prop) %>% 
  summarize(gene_count = sum(mean_bin > MEAN_BIN_THRESHOLD), 
            sig_gene_count = sum(mean_bin > MEAN_BIN_THRESHOLD &  fdr.x < FDR_THRESHOLD), 
            tasic_total_count = sum(mean_bin > MEAN_BIN_THRESHOLD & fdr < FDR_THRESHOLD), 
            tasic_overlap_count = sum(mean_bin > MEAN_BIN_THRESHOLD & fdr.x < FDR_THRESHOLD & fdr < FDR_THRESHOLD & cons_weighted == T), 
            tasic_overlap_pct = tasic_overlap_count / sig_gene_count, 
            overlap_p = 1 - phyper(tasic_overlap_count, tasic_total_count, gene_count - tasic_total_count, sig_gene_count))%>% 
  mutate(weighted_val = "corrected")

patch_seq_tasic_consistency_unweighted = merged_ephys_gene_corrs_w_aibs %>% group_by(dataset, ephys_prop) %>% 
  summarize(gene_count = sum(mean_bin > MEAN_BIN_THRESHOLD), 
            sig_gene_count = sum(mean_bin > MEAN_BIN_THRESHOLD &  fdr.y < FDR_THRESHOLD), 
            tasic_total_count = sum(mean_bin > MEAN_BIN_THRESHOLD & fdr < FDR_THRESHOLD), 
            tasic_overlap_count = sum(mean_bin > MEAN_BIN_THRESHOLD & fdr.y < FDR_THRESHOLD & fdr < FDR_THRESHOLD & cons_unweighted == T), 
            tasic_overlap_pct = tasic_overlap_count / sig_gene_count, 
            overlap_p = 1 - phyper(tasic_overlap_count, tasic_total_count, gene_count - tasic_total_count, sig_gene_count)) %>% 
  mutate(weighted_val = "uncorrected")

patch_seq_tasic_consistency = bind_rows(patch_seq_tasic_consistency, patch_seq_tasic_consistency_unweighted)


patch_seq_tasic_consistency$plabel = paste('p =\n', formatC(patch_seq_tasic_consistency$overlap_p, format = "f", digits = 3))
# patch_seq_tasic_consistency$plabel_unweighted = paste('p =\n', formatC(patch_seq_tasic_consistency$overlap_p_unweighted, format = "f", digits = 3))
# patch_seq_tasic_consistency[patch_seq_tasic_consistency$overlap_p < .01, 'plabel'] = formatC(patch_seq_tasic_consistency[patch_seq_tasic_consistency$overlap_p < .01, 'overlap_p'], format = "e", digits = 2)

patch_seq_tasic_consistency_filt = patch_seq_tasic_consistency %>% filter(tasic_total_count > 0)
# patch_seq_tasic_consistency_filt = patch_seq_tasic_consistency %>% filter(sig_gene_count > 0, tasic_total_count > 0)

patch_seq_tasic_consistency_filt$ephys_prop = plyr::mapvalues(patch_seq_tasic_consistency_filt$ephys_prop, from = response_vars, to = response_vars_plot_names)

dataset_formal_names = c('Cadwell', 'Földy', 'Fuzik', 'Tasic/AIBS')
dataset_unformal_names = c('cadwell', 'foldy', 'fuzik', 'aibs')
patch_seq_tasic_consistency_filt$dataset = plyr::mapvalues(patch_seq_tasic_consistency_filt$dataset, from = dataset_unformal_names, to = dataset_formal_names)
patch_seq_tasic_consistency_filt$weighted_val = factor(patch_seq_tasic_consistency_filt$weighted_val, levels= c('uncorrected', 'corrected'))
show_ephys_vars = c('ahpamp', 'rmp', 'aphw')
show_ephys_vars = plyr::mapvalues(show_ephys_vars, from = response_vars, to = response_vars_plot_names)

patch_seq_aibs_consistency_fig = patch_seq_tasic_consistency_filt %>% filter(ephys_prop %in% show_ephys_vars) %>% 
  ggplot(aes(x = dataset, y = tasic_overlap_count, label = plabel, fill = weighted_val)) + 
  geom_bar(stat = 'identity', position = "dodge")+ facet_wrap(~ephys_prop) + 
  # geom_text(size = 4, position = position_dodge(width = 1), color = "red") +
  ylab('Sig. gene set overlap \nwith AIBS/Tasic (count)') + xlab('Patch-seq dataset')+ 
  scale_fill_manual("", values = c("uncorrected" = "grey50", "corrected" = "black"))
patch_seq_aibs_consistency_fig
ggsave(filename = 'plots/patch_seq_aibs_consistency_fig.pdf', plot = patch_seq_aibs_consistency_fig, units = "in", height = 4, width = 10)


### put together basic example figure illustrating correlations across datasets and in aibs
gene = 'Nek7'
ephys_prop = 'aphw'

merged_ephys_gene_corrs_w_aibs %>% filter(gene == 'Nek7', ephys_prop == 'aphw')

joined_df = patch_seq_datasets$foldy$joined_df

use_cell_ids = joined_df[, c('major_type', 'cell_id', 'num_genes', 'contam_sum', 'marker_sum_norm', 'rmp', 'aphw')] %>% select(cell_id) %>% unlist

use_dataset_inds = which(joined_df$cell_id %in% use_cell_ids & !is.na(joined_df$aphw))
# weights = 1/(lapply(joined_df$contam_sum[use_dataset_inds], function(x) max(x, .2)) %>% unlist)
weights = lapply(joined_df$dissoc_corr[use_dataset_inds], function(x) max(x, MIN_QUALITY_VALUE)) %>% unlist

joined_df = patch_seq_datasets$foldy$joined_df[use_dataset_inds, c(gene, ephys_prop, 'dissoc_corr', 'colors')]
joined_df$dataset = 'Földy (corrected)'
joined_df$weights = weights %>% scale(., center = F, scale = T)
joined_df_total = joined_df

joined_df = patch_seq_datasets$foldy$joined_df[use_dataset_inds, c(gene, ephys_prop, 'dissoc_corr', 'colors')]
joined_df$weights = rep.int(1, length(weights)) / mean(weights) + runif(length(weights), 0, .0001)
joined_df$weights = joined_df$weights  %>% scale(., center = F, scale = T)
joined_df$dataset = 'Földy (uncorrected)'
joined_df_total = bind_rows(joined_df_total, joined_df)


ephys_label = getEphysNiceName(ephys_prop)


joined_df = aibs_cre_joined[ , c(gene, ephys_prop, 'sample_size_weights', 'colors')]
joined_df$dataset = 'AIBS/Tasic (pooled)'
joined_df$weights = joined_df$sample_size_weights / mean(joined_df$sample_size_weights )
joined_df$weights = joined_df$weights  %>% scale(., center = F, scale = T)
joined_df_total = bind_rows(joined_df_total, joined_df)

joined_df_total$colors = 'black'

ephys_label = getEphysNiceName(ephys_prop)

joined_df_total$dataset = factor(joined_df_total$dataset, c('Földy (uncorrected)', 'Földy (corrected)', 'AIBS/Tasic (pooled)'))

ephys_expr_example_fig = ggplot(joined_df_total, 
            aes_string(x = gene, y = ephys_prop, size = 'weights', color = 'colors')) +
  # geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') +
  geom_smooth(method = "lm", aes(weight = weights), se = F, linetype = 1, color = 'black') +
  geom_point(alpha = .5)  + 
  # coord_cartesian(xlim = c(min(joined_df[gene])-1, max(joined_df[gene])), ylim, expand = F) + 
  scale_x_continuous(trans = 'log2') + annotation_logticks(sides = "bl") + 
  scale_color_identity() + 
  scale_y_log10(breaks = c(.5, 1, 2)) + 
  theme(legend.position="none")  + facet_wrap(~dataset)
ephys_expr_example_fig = ephys_expr_example_fig + xlab(paste(gene, '(TPM + 1)') )  + ylab(ephys_label)
ephys_expr_example_fig
corr_fig = plot_grid(ephys_expr_example_fig, patch_seq_gene_count_fig, nrow = 1)

# make big combined figure


combined_fig_2 = plot_grid(gene_count_fig, corr_fig, nrow = 2, rel_heights = c(.8, 1))
ggsave(filename = 'plots/comb_gene_count_ephys_fig_v2.pdf', plot = combined_fig_2, units = "in", height = 7, width = 13)



