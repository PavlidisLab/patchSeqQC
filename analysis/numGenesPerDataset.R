library(scales)
library(ggrepel)
print('Counting numbers of genes expressed per patch seq dataset')

### set up num genes per dataset figure

select_fields = c('cell_id', 'num_genes', 'major_type', 'colors', 'ercc_pct', 'contam_sum', 'marker_sum_norm', 'read_count', 'rmp')

dfs_meta = lapply(names(patch_seq_datasets), function(dataset){
  df = patch_seq_datasets[[dataset]]$joined_df %>% select(one_of(select_fields)) %>% mutate(dataset = dataset)
  return(df)
}) %>% bind_rows()

dfs_meta$dataset = mapvalues(dfs_meta$dataset, from = c("cadwell", "foldy", "fuzik"), to = c("Cadwell", "Földy", "Fuzik"))

## generate figure showing num genes per dataset and spike in ratio and contamination for cadwell dataset

p2 = ggplot(data = dfs_meta, aes(x = dataset, y = read_count)) + geom_violin() + scale_y_log10(labels = comma)
# plot_grid(p1, p2)

joined_df = patch_seq_datasets$cadwell$joined_df
p0 = ggplot(joined_df, aes(x = ercc_pct, y = num_genes)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = .6) + 
  scale_color_identity()  +  
  xlab('Spike-in mRNA ratio (%)') + ylab('Detected genes')
p1 = ggplot(data = dfs_meta, aes(x = dataset, y = num_genes)) + geom_violin() + ylab('Detected genes') + xlab('Dataset')
p2 = ggplot(joined_df, aes(x = contam_sum, y = num_genes)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = .6) + scale_color_identity() + 
  xlab('Contamination score (AU)') + ylab('Detected genes')


## calculate ercc percents for tasic data
tasic_df = aibsExprDataDf %>% select(one_of(c('sample_title', 'read_count', 'num_genes', 'ercc_counts', 'td_tomato_pct')))
rownames(tasic_df) = tasic_df$sample_title
tasic_df = tasic_df %>% mutate(ercc_pct = 100 * ercc_counts/read_count) 

tasic_df = cbind(tasic_df, aibs_contam_all_sub)
tasic_ndnf_num_genes_df = tasic_df %>% filter(contam_type == 'Ndnf_on')

levels(aibs_med_exprs$contam_type) <- c(levels(aibs_med_exprs$contam_type), "Ndnf_on")
aibs_med_exprs$contam_type[aibs_med_exprs$contam_type == 'Ndnf'] <- 'Ndnf_on'
aibs_med_exprs$contam_type
contam_sum = normalizeContamSum(tasic_ndnf_num_genes_df, aibs_med_exprs, compare_cell_types_inh)
tasic_ndnf_num_genes_df$contam_sum = contam_sum

cadwell_df = patch_seq_datasets$cadwell$joined_df[patch_seq_datasets$cadwell$joined_df$major_type %in% c('eNGC', 'SBC'), c('num_genes', 'read_count', 'ercc_pct', 'contam_sum')]
cadwell_model = lm(num_genes ~ log10(read_count) + log10(ercc_pct) +contam_sum, data = cadwell_df)
tasic_ndnf_model = lm(num_genes ~ log10(read_count) + log10(ercc_pct) +contam_sum, data = tasic_ndnf_num_genes_df)

afss <- anova(tasic_ndnf_model)$"Sum Sq"
rsq = summary(tasic_ndnf_model)$r.squared

var_exp_df_tasic = data.frame(variable = c('Library size', 'Spike-in ratio', 'Contam score', 'Residual'))
var_exp_df_tasic = cbind(var_exp_df_tasic, percent_exp =(afss/sum(afss))/rsq*100, var_exp = (afss/sum(afss)) * 100)
var_exp_df_tasic$dataset = 'Tasic Ndnf'

afss <- anova(cadwell_model)$"Sum Sq"
rsq = summary(cadwell_model)$r.squared

var_exp_df_cadwell = data.frame(variable = c('Library size', 'Spike-in ratio', 'Contam score', 'Residual'))
var_exp_df_cadwell = cbind(var_exp_df_cadwell, percent_exp =(afss/sum(afss))/rsq*100, var_exp = (afss/sum(afss)) * 100)
var_exp_df_cadwell$dataset = 'Cadwell'

var_exp_comb = rbind(var_exp_df_tasic, var_exp_df_cadwell)
var_exp_comb$dataset = factor(var_exp_comb$dataset, levels = c('Tasic Ndnf', 'Cadwell'))
# var_exp_comb$dataset = factor(var_exp_comb$dataset, levels = c('Tasic Ndnf', 'Cadwell'))

# drop residual term
var_exp_comb = var_exp_comb[var_exp_comb$variable != 'Residual', ]

var_exp_comb$variable = factor(var_exp_comb$variable, levels = c('Library size', 'Spike-in ratio', 'Contam score'))

# make a figure that summarizes variance explained
var_exp_fig = ggplot(var_exp_comb, aes(x = variable, y = var_exp, fill = dataset)) + geom_bar(stat = 'identity', position = 'dodge') + 
  ylab('Var. Explained (norm. %)')+ theme(legend.position="top") + xlab('') + 
  scale_fill_manual("", values = c("Tasic Ndnf" = "grey80", "Cadwell" = "grey30"))

# create a combined figure that combines all variance explained fig subparts
gene_count_fig = plot_grid(p1, p0, p2, var_exp_fig, nrow = 1)
gene_count_fig


### plot spike in ratio and cell contamination measures for foldy dataset

joined_df = patch_seq_datasets$foldy$joined_df[patch_seq_datasets$foldy$joined_df$ercc_pct > 0, ]
show_cells = c('FS.INT.3', 'FS.INT.7', 'FS.INT.2')

joined_df$cell_id[!joined_df$cell_id %in% show_cells] = '' 
p0 = ggplot(joined_df, aes(x = ercc_pct, y = num_genes, label = cell_id)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = .6) + 
  geom_text_repel() + 
  scale_color_identity()  +  
  xlab('Spike-in mRNA ratio (%)') + ylab('Detected genes')
p2 = ggplot(joined_df, aes(x = contam_sum, y = num_genes, label = cell_id)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = .6) + 
  geom_text_repel() + 
  scale_color_identity() + 
  xlab('Cell contamination (AU)') + ylab('Detected genes') 
p3 = ggplot(joined_df, aes(x = ercc_pct, y = contam_sum, label = cell_id)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = .6) + 
  geom_text_repel() + 
  scale_color_identity() + 
  xlab('Spike-in mRNA ratio (%)') + ylab('Cell contamination (AU)') 
p4 = ggplot(joined_df, aes(x = ercc_pct, y = marker_sum_norm, label = cell_id)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = .6) + 
  geom_text_repel() + 
  scale_color_identity() + 
  xlab('Spike-in mRNA ratio (%)') + ylab('Endogenous marker ratio') 

foldy_gene_count_fig = plot_grid(p0, p2, p3, p4, nrow = 2, align = 'hv')
foldy_gene_count_fig

patch_seq_datasets$foldy$joined_df[, 
                                   c('cell_id', 'ercc_pct', 'num_genes', 'dissoc_corr', 'marker_sum_norm', 'major_type', 'contam_sum')] %>% 
  filter(ercc_pct > 0) 


cell_type_colors = c(Ndnf_on = "indianred1", Pyramidal = "turquoise", Pyramidal_on = "turquoise", 
                     Microglia = "gray50", Sncg_on = "purple", Astrocyte = "green4", Pvalb = "red", Pvalb_on = "red",
                     Inhibitory = "red",
                     Endothelial = "brown", 
                     Oligodendrocyte = "sandybrown",
                     OPC = 'darkorange4') %>% as.data.frame %>% tibble::rownames_to_column(var = "cell_type")
colnames(cell_type_colors) = c('cell_type', 'colors')

contam_scores = patch_seq_datasets$foldy$contam_scores
contam_scores[contam_scores$cell_id %in% show_cells, 'dissoc_corr']



fullMarkerListPvalb = fullMarkerList[c('Pvalb_on', compare_cell_types_inh)]
marker_df_pvalb = do.call(rbind, lapply(seq_along(fullMarkerListPvalb), function(i){
  data.frame(CLUSTER=names(fullMarkerListPvalb)[i], fullMarkerListPvalb[[i]])
}))

colnames(marker_df_pvalb) = c('cell_type', 'gene')


markers_collapsed = fullMarkerList[c('Pvalb_on', compare_cell_types_inh)] %>% unlist %>% unique
aibs_pvalb = aibsExprDataDf[, c('contam_type', markers_collapsed)] %>% filter(contam_type == 'Pvalb_on')

aibs_pvalb_mean = aibs_pvalb %>% select(-contam_type) %>% log2 %>% colMeans
aibs_pvalb_mean = 2^aibs_pvalb_mean

foldy_expr = patch_seq_datasets$foldy$joined_df[, c(markers_collapsed)] 
rownames(foldy_expr) = patch_seq_datasets$foldy$joined_df$cell_id

foldy_aibs_df = cbind(t(foldy_expr), aibs_pvalb_mean) %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% left_join(., marker_df_pvalb, by = 'gene') %>% left_join(., cell_type_colors)

foldy_aibs_df = foldy_aibs_df %>% gather(cell_name, expr_val, FS.INT:FS.INT.9)

# show_cells = c('RS.INT.12', 'RS.INT.11', 'RS.INT.10')
foldy_aibs_df_plot= foldy_aibs_df %>% filter(cell_name %in% show_cells)
foldy_aibs_df_plot$cell_name = factor(foldy_aibs_df_plot$cell_name, levels = show_cells)

pvalb_marker_illustration_plot = foldy_aibs_df_plot  %>% ggplot(aes(x = expr_val, aibs_pvalb_mean, color = colors)) + 
  geom_smooth(method = "lm", se = F, linetype = 2, color = 'black') + 
  geom_point() + scale_y_continuous(trans = "log2") + scale_x_continuous(trans = "log2") + geom_abline(slope = 1, intercept = 0) +
  scale_color_identity("Cell type markers", labels = foldy_aibs_df_plot$cell_type, breaks = foldy_aibs_df_plot$colors, guide = "legend") +
  xlab('Patch-seq expr (TPM + 1)') + ylab('Dissoc. cell expr (TPM + 1)') + 
  facet_wrap(~cell_name)

foldy_comb_gene_count_fig = plot_grid(foldy_gene_count_fig,pvalb_marker_illustration_plot, nrow =2, rel_heights  = c(2, 1))

ggsave(filename = 'plots/foldy_comb_gene_count_fig.pdf', plot = foldy_comb_gene_count_fig, units = "in", height = 10, width = 8)

### plot contamination vs normalized marker sum cell contamination measures for cadwell and foldy dataset

contam_scores = patch_seq_datasets$cadwell$contam_scores %>% filter(major_type %in% c('eNGC', 'SBC'))
eval_cor = cor(contam_scores$marker_sum, contam_scores$contam_sum)

show_cell_labels_cadwell = c('eNGC.13', 'eNGC.12')
contam_scores$cell_id[!contam_scores$cell_id %in% show_cell_labels_cadwell] = ''
contam_scores[contam_scores$cell_id %in% show_cell_labels_cadwell, 'dissoc_corr']

p0 = ggplot(contam_scores, aes(x = marker_sum, y = contam_sum, label = cell_id)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = 1) + xlab('Ndnf markers (log2 TPM+1)') + ylab('Cell contamination (AU)') + 
  geom_text_repel(fontface = 'bold') + 
  ggtitle(paste('Cadwell Ndnf cells,\nr = ', sprintf("%.2f", eval_cor)))

contam_scores = patch_seq_datasets$foldy$contam_scores %>% filter(major_type %in% c('RS-INT'))
eval_cor = cor(contam_scores$marker_sum, contam_scores$contam_sum)

show_cell_labels_foldy = c('RS.INT.12', 'RS.INT.3')
contam_scores[contam_scores$cell_id %in% show_cell_labels_foldy, 'dissoc_corr']
contam_scores$cell_id[!contam_scores$cell_id %in% show_cell_labels_foldy] = ''

p2 = ggplot(contam_scores, aes(x = marker_sum, y = contam_sum, label = cell_id)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
  geom_point(size = 2, alpha = 1) + xlab('Sncg markers (log2 TPM+1)') + ylab('Cell contamination (AU)') +
  geom_text_repel(fontface = 'bold') + 
  ggtitle(paste('Földy Sncg cells,\nr = ', sprintf("%.2f", eval_cor)))
# p2 = ggplot(joined_df, aes(x = contam_sum, y = num_genes)) + geom_smooth(method = "lm", se = F, linetype = 1, color = 'grey') + 
#   geom_point(size = 2, alpha = .6) + scale_color_identity() + 
#   xlab('Cell contamination (AU)') + ylab('Detected genes') + ylim(0, 12500)


### illustrate calculation of quality scores using correlation with dissociated cell data along markers
fullMarkerListNdnf = fullMarkerList[c('Ndnf_on', compare_cell_types_inh)]
marker_df_ndnf = do.call(rbind, lapply(seq_along(fullMarkerListNdnf), function(i){
  data.frame(CLUSTER=names(fullMarkerListNdnf)[i], fullMarkerListNdnf[[i]])
}))
colnames(marker_df_ndnf) = c('cell_type', 'gene')

markers_collapsed = fullMarkerList[c('Ndnf_on', compare_cell_types_inh)] %>% unlist %>% unique
aibs_ndnf = aibsExprDataDf[, c('contam_type', markers_collapsed)] %>% filter(contam_type == 'Ndnf_on')

aibs_ndnf_mean = aibs_ndnf %>% select(-contam_type) %>% log2 %>% colMeans
aibs_ndnf_mean = 2^aibs_ndnf_mean

cadwell_expr = patch_seq_datasets$cadwell$joined_df[, c(markers_collapsed)] 
rownames(cadwell_expr) = patch_seq_datasets$cadwell$joined_df$cell_id

cadwell_aibs_df = cbind(t(cadwell_expr), aibs_ndnf_mean) %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  left_join(., marker_df_ndnf, by = 'gene') %>% left_join(., cell_type_colors)

cadwell_aibs_df = cadwell_aibs_df %>% gather(cell_name, expr_val, eNGC:eNGC.22)

cadwell_aibs_df_plot= cadwell_aibs_df %>% filter(cell_name %in% show_cell_labels_cadwell)
cadwell_aibs_df_plot$cell_name = factor(cadwell_aibs_df_plot$cell_name, levels = show_cell_labels_cadwell)


cadwell_marker_illust_plot = cadwell_aibs_df_plot  %>% ggplot(aes(x = expr_val, aibs_ndnf_mean, color = colors)) + 
  geom_smooth(method = "lm", se = F, linetype = 2, color = 'black') + 
  geom_point() + scale_y_continuous(trans = "log2") + scale_x_continuous(trans = "log2") + geom_abline(slope = 1, intercept = 0) +
  scale_color_identity("Cell type markers", labels = cadwell_aibs_df_plot$cell_type, breaks = cadwell_aibs_df_plot$colors, guide = "legend") +
  xlab('Patch-seq expr (TPM + 1)') + ylab('Dissoc. cell expr (TPM + 1)') + 
  facet_wrap(~cell_name)

fullMarkerListSncg = fullMarkerList[c('Sncg_on', compare_cell_types_inh)]
marker_df_sncg = do.call(rbind, lapply(seq_along(fullMarkerListSncg), function(i){
  data.frame(CLUSTER=names(fullMarkerListSncg)[i], fullMarkerListSncg[[i]])
}))
colnames(marker_df_sncg) = c('cell_type', 'gene')


markers_collapsed = fullMarkerList[c('Sncg_on', compare_cell_types_inh)] %>% unlist %>% unique
aibs_sncg = aibsExprDataDf[, c('contam_type', markers_collapsed)] %>% filter(contam_type == 'Sncg_on')

aibs_sncg_mean = aibs_sncg %>% select(-contam_type) %>% log2 %>% colMeans
aibs_sncg_mean = 2^aibs_sncg_mean

foldy_expr = patch_seq_datasets$foldy$joined_df[, c(markers_collapsed)] 
rownames(foldy_expr) = patch_seq_datasets$foldy$joined_df$cell_id

foldy_aibs_df = cbind(t(foldy_expr), aibs_sncg_mean) %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% left_join(., marker_df_sncg, by = 'gene') %>% left_join(., cell_type_colors)

foldy_aibs_df = foldy_aibs_df %>% gather(cell_name, expr_val, RS.INT:RS.INT.9)

foldy_aibs_df_plot= foldy_aibs_df %>% filter(cell_name %in% show_cell_labels_foldy)
foldy_aibs_df_plot$cell_name = factor(foldy_aibs_df_plot$cell_name, levels = show_cell_labels_foldy)

foldy_marker_illust_plot = foldy_aibs_df_plot  %>% ggplot(aes(x = expr_val, aibs_sncg_mean, color = colors)) + 
  geom_smooth(method = "lm", se = F, linetype = 2, color = 'black') + 
  geom_point() + scale_y_continuous(trans = "log2") + scale_x_continuous(trans = "log2") + geom_abline(slope = 1, intercept = 0) +
  scale_color_identity("Cell type markers", labels = foldy_aibs_df_plot$cell_type, breaks = foldy_aibs_df_plot$colors, guide = "legend") +
  xlab('Patch-seq expr (TPM + 1)') + ylab('Dissoc. cell expr (TPM + 1)') + 
  facet_wrap(~cell_name)


cadwell_supp = plot_grid(p0, cadwell_marker_illust_plot, nrow = 1, rel_widths = c(1, 2.5))
foldy_supp = plot_grid(p2, foldy_marker_illust_plot, nrow = 1, rel_widths = c(1, 2.5))
cell_contam_vs_marker_expr = plot_grid(cadwell_supp, foldy_supp, nrow = 2)

# cell_contam_vs_marker_expr = plot_grid(p0, p2, nrow = 1)
ggsave(filename = 'plots/cell_contam_vs_marker_expr.pdf', plot = cell_contam_vs_marker_expr, units = "in", height = 8, width = 12)


