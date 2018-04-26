## illustrates issue of off cell type contamination in patch seq datasets

# show marker gene contamination as heatmap

source('R/calcCellContam.R')
source('R/plotMarkerHeatmaps.R')


cadwell_expr = patch_seq_datasets$cadwell$joined_df[patch_seq_datasets$cadwell$joined_df$major_type == 'eNGC', ] 
rownames(cadwell_expr)  = patch_seq_datasets$cadwell$joined_df[patch_seq_datasets$cadwell$joined_df$major_type == 'eNGC', 'cell_id']

marker_sum = patch_seq_datasets$cadwell$joined_df[patch_seq_datasets$cadwell$joined_df$major_type == 'eNGC', ]$dissoc_corr
# on_marker_sum = sumExpression(cadwell_expr, ndnf_markers)
sorted_cell_ids = order(-marker_sum)
cadwell_expr =cadwell_expr[sorted_cell_ids, ] 

plot_cell_types_cadwell = c('Ndnf_on', 'Pyramidal', 'Microglia')

plot_marker_list = c(onCellMarkerGenes[plot_cell_types_cadwell[1]], offCellMarkerGenes[plot_cell_types_cadwell[2:3]])

cadwell_contam_heatmap = plotMarkerHeatmap(plot_marker_list, cadwell_expr, show_legend = F)

## plot sum of marker expression for patch seq relative to aibs dissociated cells

cadwell_contam = calcContamAllTypes(cadwell_expr, fullMarkerList) 
cadwell_contam$colors = "indianred1"
cadwell_contam$major_type = cadwell_expr$major_type
cadwell_contam$contam_type = "Ndnf_on"
# cadwell_contam = left_join(cadwell_contam, scde_list$cadwell$qc_df %>% select(cell_id, mt_pct, num_genes))


aibs_ndnf_contam = calcContamAllTypes(aibsExprDataDf[str_detect(aibsExprDataDf$primary_type, 'Ndnf'), ], fullMarkerList) 
aibs_ndnf_contam$colors = 'indianred1'

aibs_vip_contam = calcContamAllTypes(aibsExprDataDf[str_detect(aibsExprDataDf$primary_type, 'Sncg'), ], fullMarkerList) 
aibs_vip_contam$colors = "purple"

aibs_pyr_contam = calcContamAllTypes(aibsExprDataDf[str_detect(aibsExprDataDf$primary_type, 'L'), ], fullMarkerList) 
aibs_pyr_contam$colors = "turquoise"

aibs_micro_contam = calcContamAllTypes(aibsExprDataDf[str_detect(aibsExprDataDf$primary_type, 'Micro'), ], fullMarkerList) 
aibs_micro_contam$colors = "gray50"
# aibs_micro_contam$Microglia[aibs_micro_contam$Microglia > 110] = 110

aibs_contam = rbind(aibs_ndnf_contam %>% select(Ndnf_on, Sncg_on, Microglia, Pyramidal, colors), 
                    aibs_pyr_contam %>% select(Ndnf_on, Sncg_on,  Microglia, Pyramidal, colors), 
                    aibs_vip_contam %>% select(Ndnf_on, Sncg_on, Microglia, Pyramidal, colors), 
                    aibs_micro_contam %>% select(Ndnf_on, Sncg_on, Microglia, Pyramidal, colors)
)

x_min = 0
x_max = max(aibs_contam$Ndnf_on, cadwell_contam$Ndnf_on) + 0
y_min_1 = min(c(aibs_contam$Pyramidal, cadwell_contam$Pyramidal, aibs_contam$Pyramidal)) - 0
y_max_1 = max(c(aibs_contam$Pyramidal, cadwell_contam$Pyramidal, aibs_contam$Pyramidal)) + 0

y_min_2 = min(c(aibs_contam$Microglia, cadwell_contam$Microglia, aibs_contam$Microglia)) - 0
y_max_2 = max(aibs_contam$Microglia, cadwell_contam$Microglia, aibs_contam$Microglia) + 0

PLOT_THRESHOLD = .025

aibs_pyr = aibs_contam %>% filter(colors == 'turquoise') %>% select(Pyramidal) %>% unlist
aibs_pyr_low = aibs_pyr %>% quantile(PLOT_THRESHOLD)
aibs_pyr_high = aibs_pyr %>% quantile(1-PLOT_THRESHOLD)
aibs_ndnf = aibs_contam %>% filter(colors == 'indianred1') %>% select(Ndnf_on) %>% unlist
aibs_ndnf_low = aibs_ndnf %>% quantile(PLOT_THRESHOLD)
aibs_ndnf_high = aibs_ndnf %>% quantile(1-PLOT_THRESHOLD)
aibs_micro = aibs_contam %>% filter(colors == 'gray50') %>% select(Microglia) %>% unlist
aibs_micro_low = aibs_micro %>% quantile(PLOT_THRESHOLD)
aibs_micro_high = aibs_micro %>% quantile(1-PLOT_THRESHOLD)
aibs_vip = aibs_contam %>% filter(colors == 'purple') %>% select(Sncg_on) %>% unlist
aibs_vip_low = aibs_vip %>% quantile(PLOT_THRESHOLD)
aibs_vip_high = aibs_vip %>% quantile(1-PLOT_THRESHOLD)

p1 = ggplot(aibs_contam %>% filter(colors %in% c('turquoise', 'indianred1')), aes(x = Ndnf_on, y = Pyramidal, color = colors)) + 
  geom_hline(yintercept = aibs_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
  geom_point(alpha = 1) + scale_color_identity() + 
  xlab('Ndnf markers (log2 TPM+1)') + ylab('Pyr markers (log2 TPM+1)') +
  ggtitle('Ndnf and Pyr cells,\nTasic (dissociated)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12)) 

  
p2 = ggplot(cadwell_contam, aes(x = Ndnf_on, y = Pyramidal, color = colors)) + 
  geom_hline(yintercept = aibs_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() + 
  xlab('Ndnf markers (log2 TPM+1)') + ylab('Pyr markers (log2 TPM+1)') + 
  ggtitle('eNGC cells,\nCadwell (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p3 = ggplot(aibs_contam %>% filter(colors %in% c('indianred1', 'gray50')), aes(x = Ndnf_on, y = Microglia, color = colors)) + 
  geom_hline(yintercept = aibs_micro_high, color = 'gray50', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_micro_low, color = 'gray50', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
  geom_point(alpha =1 ) + scale_color_identity() + 
  xlab('Ndnf markers (log2 TPM+1)') + ylab('Microglia markers (log2 TPM+1)') +
  ggtitle('Ndnf and Microglia cells,\nTasic (dissociated)') + xlim(x_min, x_max) + ylim(y_min_2, y_max_2) + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))
p4 = ggplot(cadwell_contam, aes(x = Ndnf_on, y = Microglia, color = colors)) + 
  geom_hline(yintercept = aibs_micro_high, color = 'gray50', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_micro_low, color = 'gray50', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() + 
  xlab('Ndnf markers (log2 TPM+1)') + ylab('Microglia markers (log2 TPM+1)') +
  ggtitle('eNGC cells,\nCadwell (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_2, y_max_2) + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

cadwell_aibs_marker_comp = plot_grid(p1, p2, p3, p4, nrow = 2)

cadwell_contam_combined = plot_grid(cadwell_contam_heatmap$gtable, cadwell_aibs_marker_comp, nrow = 1, rel_widths = c(4, 6))


### make similar plots for foldy RS cells data


# show marker gene contamination as heatmap

foldy_expr = patch_seq_datasets$foldy$joined_df[patch_seq_datasets$foldy$joined_df$major_type == 'RS-INT', ] 
rownames(foldy_expr)  = foldy_expr$cell_id

marker_sum = patch_seq_datasets$foldy$joined_df[patch_seq_datasets$foldy$joined_df$major_type == 'RS-INT', ]$dissoc_corr
# on_marker_sum = sumExpression(cadwell_expr, ndnf_markers)
sorted_cell_ids = order(-marker_sum)
foldy_expr =foldy_expr[sorted_cell_ids, ] 

plot_cell_types_foldy = c('Sncg_on', 'Pyramidal', 'Microglia')

plot_marker_list = c(onCellMarkerGenes[plot_cell_types_foldy[1]], offCellMarkerGenes[plot_cell_types_foldy[2:3]])

foldy_contam_heatmap = plotMarkerHeatmap(plot_marker_list, foldy_expr, show_legend = T, show_cell_labels = F)

## plot sum of marker expression for patch seq relative to aibs dissociated cells

foldy_contam_vip = calcContamAllTypes(foldy_expr, fullMarkerList) 
foldy_contam_vip$colors = "purple"
foldy_contam_vip$major_type = foldy_expr$major_type
foldy_contam_vip$contam_type = 'Sncg_on'
# cadwell_contam = left_join(cadwell_contam, scde_list$cadwell$qc_df %>% select(cell_id, mt_pct, num_genes))

x_min = 0
x_max = max(aibs_vip_contam$Sncg, foldy_contam_vip$Sncg_on) + 0
y_min_1 = min(c(aibs_contam$Pyramidal, foldy_contam_vip$Pyramidal)) - 0
y_max_1 = max(c(aibs_contam$Pyramidal, foldy_contam_vip$Pyramidal)) + 0

y_min_2 = min(c(aibs_contam$Microglia, foldy_contam_vip$Microglia)) - 0
y_max_2 = max(c(aibs_contam$Microglia, foldy_contam_vip$Microglia)) + 0

p1 = ggplot(aibs_contam %>% filter(colors %in% c('turquoise', 'purple')), aes(x = Sncg_on, y = Pyramidal, color = colors)) + 
  geom_hline(yintercept = aibs_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_vip_high, color = 'purple', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_vip_low, color = 'purple', linetype = 2, alpha = .5) +
  geom_point(alpha = 1) + scale_color_identity() + 
  xlab('Sncg markers (log2 TPM+1)') + ylab('Pyr markers (log2 TPM+1)') + 
  ggtitle('Sncg and Pyr cells,\nTasic (dissociated)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p3 = ggplot(foldy_contam_vip, aes(x = Sncg_on, y = Pyramidal, color = colors)) + 
  geom_hline(yintercept = aibs_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_vip_high, color = 'purple', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_vip_low, color = 'purple', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() + 
  xlab('Sncg markers (log2 TPM+1)') + ylab('Pyr markers (log2 TPM+1)') + 
  ggtitle('RS-INT cells, \nFöldy (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p4 = ggplot(aibs_contam %>% filter(colors %in% c('purple', 'gray50')), aes(x = Sncg_on, y = Microglia, color = colors)) + 
  geom_hline(yintercept = aibs_micro_high, color = 'gray50', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_micro_low, color = 'gray50', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_vip_high, color = 'purple', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_vip_low, color = 'purple', linetype = 2, alpha = .5) +
  geom_point(alpha = 1) + scale_color_identity() + 
  xlab('Sncg markers (log2 TPM+1)') + ylab('Microglia markers (log2 TPM+1)') + 
  ggtitle('Sncg and Microglia cells,\nTasic (dissociated)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_2) +       
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p6 = ggplot(foldy_contam_vip, aes(x = Sncg_on, y = Microglia, color = colors)) + 
  geom_hline(yintercept = aibs_micro_high, color = 'gray50', linetype = 2, alpha = .5) + geom_hline(yintercept = aibs_micro_low, color = 'gray50', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = aibs_vip_high, color = 'purple', linetype = 2, alpha = .5) + geom_vline(xintercept = aibs_vip_low, color = 'purple', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() + 
  xlab('Sncg markers (log2 TPM+1)') + ylab('Microglia markers (log2 TPM+1)') +
  ggtitle('RS-INT cells, \nFöldy (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_2) +       
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

foldy_aibs_marker_comp = plot_grid(p1, p3, p4, p6, nrow = 2)

foldy_contam_combined = plot_grid(foldy_contam_heatmap$gtable, foldy_aibs_marker_comp, nrow = 1, rel_widths = c(4, 6))


combined_contam_plot = plot_grid(cadwell_contam_combined, foldy_contam_combined, nrow = 2)

ggsave(filename = 'plots/contam_main_fig.pdf', plot = combined_contam_plot, units = "in", width = 8, height = 10)


### show heatmap for zeisel and fuzik


fuzik_expr = patch_seq_datasets$fuzik$joined_df[patch_seq_datasets$fuzik$joined_df$major_type == 'Exc', ] 
rownames(fuzik_expr)  = patch_seq_datasets$fuzik$joined_df[patch_seq_datasets$fuzik$joined_df$major_type == 'Exc', 'cell_id']

marker_sum = patch_seq_datasets$fuzik$joined_df[patch_seq_datasets$fuzik$joined_df$major_type == 'Exc', ]$dissoc_corr
# on_marker_sum = sumExpression(fuzik_expr, ndnf_markers)
sorted_cell_ids = order(-marker_sum)
fuzik_expr =fuzik_expr[sorted_cell_ids, ] 

plot_cell_types_fuzik = c('Pyramidal_on', 'Astrocyte', 'Microglia')
plot_marker_list = c(onCellMarkerGenes[plot_cell_types_fuzik[1]], offCellMarkerGenes[plot_cell_types_fuzik[2:3]])

fuzik_contam_heatmap = plotMarkerHeatmap(plot_marker_list, fuzik_expr, show_legend = T)

## plot sum of marker expression for patch seq relative to aibs dissociated cells

fuzik_contam = calcContamAllTypes(fuzik_expr, fullMarkerList) 
fuzik_contam$colors = "turquoise"
fuzik_contam$major_type = fuzik_expr$major_type
fuzik_contam$contam_type = "Pyramidal_on"
# fuzik_contam = left_join(fuzik_contam, scde_list$fuzik$qc_df %>% select(cell_id, mt_pct, num_genes))

# look at zeisel expression for just cortex data
zeisel_contam_all_sub_cortex = cbind(zeisel_contam_all_sub, zeiselExprDataDf[, c('tissue', 'level1class')]) %>% filter(tissue== 'sscortex')
zeisel_ndnf_contam =  zeisel_contam_all_sub_cortex %>% filter(contam_type == c('Ndnf_on'))
zeisel_ndnf_contam$colors = 'indianred1'

zeisel_astro_contam = zeisel_contam_all_sub_cortex %>% filter(contam_type == c('Astrocyte_on'))
zeisel_astro_contam$colors = 'green4'

zeisel_pyr_contam = zeisel_contam_all_sub_cortex %>% filter(contam_type == c('Pyramidal_on'))
zeisel_pyr_contam$colors = "turquoise"

zeisel_micro_contam = zeisel_contam_all_sub_cortex %>% filter(contam_type == c('Microglia_on'))
zeisel_micro_contam$colors = "gray50"
# zeisel_micro_contam$Microglia[zeisel_micro_contam$Microglia > 110] = 110

zeisel_contam = rbind(zeisel_ndnf_contam %>% select(Ndnf_on, Astrocyte, Microglia, Pyramidal_on, colors), 
                    zeisel_pyr_contam %>% select(Ndnf_on, Astrocyte,  Microglia, Pyramidal_on, colors), 
                    zeisel_micro_contam %>% select(Ndnf_on, Astrocyte, Microglia, Pyramidal_on, colors), 
                    zeisel_astro_contam %>% select(Ndnf_on, Astrocyte, Microglia, Pyramidal_on, colors)
)

x_min = 0
x_max = max(zeisel_contam$Pyramidal_on, fuzik_contam$Pyramidal_on) + 0
y_min_1 = min(c(zeisel_contam$Astrocyte, fuzik_contam$Astrocyte, zeisel_contam$Astrocyte)) - 0
y_max_1 = max(c(zeisel_contam$Astrocyte, fuzik_contam$Astrocyte, zeisel_contam$Astrocyte)) + 0

y_min_2 = min(c(zeisel_contam$Microglia, fuzik_contam$Microglia, zeisel_contam$Microglia)) - 0
y_max_2 = max(zeisel_contam$Microglia, fuzik_contam$Microglia, zeisel_contam$Microglia) + 0

PLOT_THRESHOLD = .025

zeisel_pyr = zeisel_contam %>% filter(colors == 'turquoise') %>% select(Pyramidal_on) %>% unlist
zeisel_pyr_low = zeisel_pyr %>% quantile(PLOT_THRESHOLD)
zeisel_pyr_high = zeisel_pyr %>% quantile(1-PLOT_THRESHOLD)
zeisel_ndnf = zeisel_contam %>% filter(colors == 'indianred1') %>% select(Ndnf_on) %>% unlist
zeisel_ndnf_low = zeisel_ndnf %>% quantile(PLOT_THRESHOLD)
zeisel_ndnf_high = zeisel_ndnf %>% quantile(1-PLOT_THRESHOLD)
zeisel_astro = zeisel_contam %>% filter(colors == 'green4') %>% select(Astrocyte) %>% unlist
zeisel_astro_low = zeisel_astro %>% quantile(PLOT_THRESHOLD)
zeisel_astro_high = zeisel_astro %>% quantile(1-PLOT_THRESHOLD)
zeisel_micro = zeisel_contam %>% filter(colors == 'gray50') %>% select(Microglia) %>% unlist
zeisel_micro_low = zeisel_micro %>% quantile(PLOT_THRESHOLD)
zeisel_micro_high = zeisel_micro %>% quantile(1-PLOT_THRESHOLD)
# zeisel_vip = zeisel_contam %>% filter(colors == 'purple') %>% select(Sncg) %>% unlist
# zeisel_vip_low = zeisel_vip %>% quantile(PLOT_THRESHOLD)
# zeisel_vip_high = zeisel_vip %>% quantile(1-PLOT_THRESHOLD)
# 
# p1 = ggplot(zeisel_contam %>% filter(colors %in% c('turquoise', 'indianred1')), aes(x = Ndnf, y = Pyramidal, color = colors)) + 
#   geom_hline(yintercept = zeisel_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = zeisel_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) + 
#   geom_vline(xintercept = zeisel_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = zeisel_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
#   geom_point(alpha = .5) + scale_color_identity() + 
#   xlab('Ndnf markers (log2 CPM+1)') + ylab('Pyr markers (log2 CPM+1)') +
#   ggtitle('Int12, Int15 and Pyr cells,\nZeisel (dissociated)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
#   theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12)) 
# 
# 
# p2 = ggplot(fuzik_contam, aes(x = Ndnf, y = Pyramidal, color = colors)) + 
#   geom_hline(yintercept = zeisel_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_hline(yintercept = zeisel_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) + 
#   geom_vline(xintercept = zeisel_ndnf_high, color = 'indianred1', linetype = 2, alpha = .5) + geom_vline(xintercept = zeisel_ndnf_low, color = 'indianred1', linetype = 2, alpha = .5) +
#   geom_point() + scale_color_identity() + 
#   xlab('Ndnf markers (log2 CPM+1)') + ylab('Pyr markers (log2 CPM+1)') + 
#   ggtitle('L1/2 interneurons,\nFuzik (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
#   theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p1 = ggplot(zeisel_contam %>% filter(colors %in% c('turquoise', 'green4')), aes(x = Pyramidal_on, y = Astrocyte, color = colors)) + 
  geom_hline(yintercept = zeisel_astro_high, color = 'green4', linetype = 2, alpha = .5) + geom_hline(yintercept = zeisel_astro_low, color = 'green4', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = zeisel_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_vline(xintercept = zeisel_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) +
  geom_point(alpha = 1) + scale_color_identity() + 
  xlab('Pyr markers (log2 CPM+1)') + ylab('Astrocyte markers (log2 CPM+1)') +
  ggtitle('Pyramidal and Astrocytes,\nZeisel (dissociated)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12)) 


p2 = ggplot(fuzik_contam, aes(x = Pyramidal_on, y = Astrocyte, color = colors)) + 
  geom_hline(yintercept = zeisel_astro_high, color = 'green4', linetype = 2, alpha = .5) + geom_hline(yintercept = zeisel_astro_low, color = 'green4', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = zeisel_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_vline(xintercept = zeisel_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() + 
  xlab('Pyr markers (log2 CPM+1)') + ylab('Astrocyte markers (log2 CPM+1)') + 
  ggtitle('Pyramidal cells,\nFuzik (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_1, y_max_1) +    
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p3 = ggplot(zeisel_contam %>% filter(colors %in% c('turquoise', 'gray50')), aes(x = Pyramidal_on, y = Microglia, color = colors)) + 
  geom_hline(yintercept = zeisel_micro_high, color = 'gray50', linetype = 2, alpha = .5) + geom_hline(yintercept = zeisel_micro_low, color = 'gray50', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = zeisel_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_vline(xintercept = zeisel_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) +
  geom_point(alpha = 1) + scale_color_identity() + 
  xlab('Pyr markers (log2 CPM+1)') + ylab('Microglia markers (log2 CPM+1)') +
  ggtitle('Pyramidal and Microglia cells,\nZeisel (dissociated)') + xlim(x_min, x_max) + ylim(y_min_2, y_max_2) + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

p4 = ggplot(fuzik_contam, aes(x = Pyramidal_on, y = Microglia, color = colors)) + 
  geom_hline(yintercept = zeisel_micro_high, color = 'gray50', linetype = 2, alpha = .5) + geom_hline(yintercept = zeisel_micro_low, color = 'gray50', linetype = 2, alpha = .5) + 
  geom_vline(xintercept = zeisel_pyr_high, color = 'turquoise', linetype = 2, alpha = .5) + geom_vline(xintercept = zeisel_pyr_low, color = 'turquoise', linetype = 2, alpha = .5) +
  geom_point() + scale_color_identity() + 
  xlab('Pyr markers (log2 CPM+1)') + ylab('Microglia markers (log2 CPM+1)') +
  ggtitle('Pyramidal cells,\nFuzik (patch-seq)') + xlim(x_min, x_max) + ylim(y_min_2, y_max_2) + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.title = element_text(size = 12))

fuzik_zeisel_marker_comp = plot_grid(p1, p2, p3, p4, nrow = 2)

fuzik_contam_combined = plot_grid(fuzik_contam_heatmap$gtable, fuzik_zeisel_marker_comp, nrow = 1, rel_widths = c(6, 6))
ggsave(filename = 'plots/fuzik_contam_combined.pdf', plot = fuzik_contam_combined, units = "in", width = 11, height = 6)

### show contamination for all single cells aggregated to the cell type level


contam_plot_types = factor(c('Astrocyte', 'Endothelial', 'Microglia', 'Oligodendrocyte', 'OPC', 'Pyramidal'), levels = c('Astrocyte', 'Endothelial', 'Microglia', 'Oligodendrocyte', 'OPC', 'Pyramidal'))
contam_mat = getContamMatrixFromPatchSeq(patch_seq_datasets$cadwell, 'eNGC', contam_show_types = contam_plot_types)
p1 = plotContamHeatmap(contam_mat)

contam_mat = getContamMatrixFromPatchSeq(patch_seq_datasets$foldy, 'RS-INT', contam_show_types = contam_plot_types)
p2 = plotContamHeatmap(contam_mat)

contam_plot_types_fuzik = factor(c('Astrocyte', 'Endothelial', 'Microglia', 'Oligodendrocyte', 'Inhibitory'), levels = c('Astrocyte', 'Endothelial', 'Microglia', 'Oligodendrocyte', 'Inhibitory'))
contam_mat = getContamMatrixFromPatchSeq(patch_seq_datasets$fuzik, 'Exc', contam_show_types = contam_plot_types_fuzik)
p3 = plotContamHeatmap(contam_mat)

all_dataset_contam = plot_grid(p1$gtable, p2$gtable, p3$gtable, nrow = 1, rel_widths = c(2, 1.8, 3))
ggsave(filename = 'plots/contam_all_cell_types.pdf', plot = all_dataset_contam, units = "in", width = 15, height = 2.5)


