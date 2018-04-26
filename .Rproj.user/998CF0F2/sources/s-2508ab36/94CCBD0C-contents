library(grid)
library(gridExtra)
library(ggrepel)

# set global variables

response_vars = c("rin", "rmp","apthr","apamp",
                  "aphw", "tau", "ahpamp", "rheo","maxfreq", "cap", "adratio")
log_vars = c("rin", "tau", "aphw", "rheo", "maxfreq", "cap")
response_vars_plot_names = c('Rin', 'Vrest', 'APthr', 'APamp', 'APhw', 'Tau', 'AHPamp', 'Rheo', 'FRmax', 'Cm', 'SFA')
response_vars_plot_names_first = c('R', 'V', 'AP', 'AP', 'AP', 'Tau', 'AHP', 'Rheo', 'FR', 'C', 'SFA')
response_vars_plot_names_second = c('in', 'rest', 'thr', 'amp', 'hw', ' ', 'amp', ' ', 'max', 'm', ' ')
response_vars_plot_names_units = c('MΩ', 'mV', 'mV', 'mV', 'ms', 'ms', 'mV', 'pA', 'Hz', 'pF', 'ratio')



ephysExprPlot = function(input_df, gene_name, ephys_name, show_labels = F, 
                         log_trans = T, aibs_flag = F, mult_plot_flag = F, remove_y_lab = F, 
                         add_trend_line = T,
                         log_trans_x = F, log_trans_x_log2 = F, 
                         plot_corr = F, corr_val, corr_pval, jitter = F, patchseq = F) {
  log_vars = c("rin", "tau", "aphw", "rheo", "maxfreq", "cap")
  
  if (mult_plot_flag){
    POINTSIZE = 2
    TITLESIZE = 10
    TEXTSIZE = 8
    LABELSIZE = 3
  } else{
    POINTSIZE = 4
    TITLESIZE = 16
    TEXTSIZE = 12
    LABELSIZE = 4
  }
  
  if (patchseq){
    POINTSIZE = 2
  }

  
  if (!(gene_name %in% colnames(input_df) & (ephys_name %in% colnames(input_df) )))
    return(plot.new())

  
  p = ggplot(input_df, aes_string(x = gene_name, y = ephys_name, color = "colors")) 
  if (add_trend_line & !patchseq){
    p = p + geom_smooth(method = "lm", aes(group=1), se = FALSE, color="grey")
  }else{
    p = p + geom_smooth(method = "lm", aes(group=1), se = FALSE, color="lightgrey", linetype = "dashed")
  }
  if (jitter== F){
    p = p + geom_point( na.rm=TRUE, size = POINTSIZE)
    p = p + scale_colour_identity() 
    p = p + geom_point(shape = 1,size = POINTSIZE,colour = "black", na.rm=TRUE) 
  }
  else{
    p = p + geom_jitter( na.rm=TRUE, size = POINTSIZE, alpha = .75, width = .05)
    p = p + scale_colour_identity() 
  }
  
  if (patchseq == T){
    mean_df1 = input_df %>% select(one_of(gene_name, 'colors')) %>% group_by(colors) %>% summarize_all(funs(2^(mean(log2(.),na.rm = TRUE))))
    mean_df2 = input_df %>% select(one_of(ephys_name, 'colors')) %>% group_by(colors) %>% summarize_all(funs(mean(.,na.rm = TRUE)))
    mean_df = cbind(mean_df1, mean_df2)
    p = p + geom_point(data = mean_df, na.rm=TRUE, size = POINTSIZE+2)
    p = p + scale_colour_identity() 
    p = p + geom_point(data = mean_df, shape = 1,size = POINTSIZE+2,colour = "black", na.rm=TRUE) 
    p = p + geom_smooth(data = mean_df, method = "lm", aes(group=2), se = FALSE, color="grey")
  }


  if (show_labels){
    if ("short.name" %in% colnames(input_df))
      p = p + geom_text_repel(aes_string(x = gene_name, y = ephys_name, label="short.name"), color = "black", na.rm = TRUE, size = LABELSIZE)
    else
      p = p + geom_text_repel(aes_string(x = gene_name, y = ephys_name, label="sample.name"), color = "black", na.rm = TRUE, size = LABELSIZE)
  }
  
  p = p + theme(axis.title = element_text(size = TITLESIZE), axis.text = element_text(size=TEXTSIZE))
  
  if (aibs_flag){
    log_trans = F
    if (ephys_name == 'maxfreq'){
      log_trans = T
    }
    gene_label = bquote(.(gene_name) ~  '(log'[2] ~'TPM+1)')
  }
  else{
    gene_label = bquote(.(gene_name) ~  '(log'[2] ~'expr)')
  }
  if ((ephys_name %in% log_vars) & (ephys_name != "aphw") & (log_trans)){
    p = p + scale_y_log10()  + annotation_logticks(sides = "l")
  }
  if ((gene_name %in% log_vars) & (gene_name != "aphw") & (log_trans)){
    p = p + scale_x_log10()  + annotation_logticks(sides = "l")
  }
  if (log_trans_x){
    p = p + scale_x_log10() + annotation_logticks(sides = "b")
  }
  if (log_trans_x_log2){
    p = p + scale_x_continuous(trans = 'log2') + annotation_logticks(sides = "b")
  }
  
  x_lower_limit = ggplot_build(p)$layout$panel_ranges[[1]]$x.range[1]

  if (!aibs_flag & x_lower_limit < 6){
    p = p + geom_vline(xintercept = 6, linetype = 2, color = 'grey', alpha = .5)
  }
  
  if (remove_y_lab){
    p = p + theme(axis.title.y = element_blank())
  }else{
    if (ephys_name %in% response_vars){
      ephys_label = getEphysNiceName(ephys_name)
      p = p + ylab(ephys_label)
    }
  }
  
  p = p + xlab(gene_label)
  
  
  if (plot_corr){
    xpos = .05
    if (corr_val < 0)
      xpos = .6
    if (aibs_flag == F){
      grob <- grobTree(textGrob(bquote(atop('r'[s] == .(corr_val),
                                            'p'[adj] == .(corr_pval))), x= xpos,  y=.85, hjust=0,
                                gp=gpar(col="black", fontsize=TEXTSIZE))
      )
    } else{
      grob <- grobTree(textGrob(bquote('r'[s] == .(corr_val)), x=xpos,  y=.95, hjust=0,
                                gp=gpar(col="black", fontsize=TEXTSIZE, fontface="italic"))
      )
    }

    p = p + annotation_custom(grob)
  }
  
  return(p)
  
}

getEphysNiceName = function(ephys_name){
  name_pos = which(response_vars == ephys_name)
  main = response_vars_plot_names_first[name_pos]
  sub = response_vars_plot_names_second[name_pos]
  unit = response_vars_plot_names_units[name_pos]
  ephys_label = paste0(response_vars_plot_names[name_pos], ' (', response_vars_plot_names_units[name_pos], ')')
  ephys_label = bquote(.(main)[.(sub)] ~  '('* .(unit) *')')
  return(ephys_label)
}


patchSeqCompPlot2 = function(gene_name, ephys_name, aibs_df = aibs_cre_joined, 
                            foldy_df = patch_seq_datasets$foldy$joined_df, 
                            cadwell_df = patch_seq_datasets$cadwell$joined_df, 
                            fuzik_df = patch_seq_datasets$fuzik$joined_df){
  # lit_df = j
  # aibs_df = jdf
  # foldy_df = foldy_j
  # cadwell_df = cadwell_joined
  
  # if (!(gene_name %in% colnames(aibs_df))){
  #   disc_gene_name = 'Hcn3'
  #   p1 = plot.new()
  # } else{
  #   disc_gene_name = gene_name
  #   p1 = ephysExprPlot(aibs_df, gene_name, ephys_name, add_trend_line=T) + ggtitle('NeuExp/NeuElec')
  # }
  
  
  # p2 = ephysExprPlot(aibs_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T)
  
  
  p2 = ggplot(data = aibs_df, aes_string(x = gene_name, y = ephys_name, size = 'sample_size_weights', weight = 'sample_size_weights', color = 'colors', group = 'broad_type')) + 
    geom_smooth(method = 'lm', aes(group=broad_type, color = broad_type_colors), se = FALSE, alpha = .5, linetype = 2) + 
    geom_smooth(method = 'lm', aes(group=1), se = FALSE, color="grey", alpha = .5, linetype = 1) + 
    geom_point() + 
    scale_color_identity() + 
    geom_point(shape = 1, na.rm=TRUE, aes(size = sample_size_weights, color = broad_type_colors), stroke = 1) +
    scale_x_continuous(trans = 'log10') + annotation_logticks(sides = "b")  + 
    theme(legend.position="none")
  p2 = p2 + xlab(paste(gene_name , '(TPM + 1)' )) + ylab(getEphysNiceName(ephys_name)) + ggtitle('AIBS')
  
  
  p3 = ephysExprPlot(cadwell_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T, jitter = T, patchseq = T)  
  if (!is.null(p3)){
    p3 = p3 + xlab(paste(gene_name , '(TPM + 1)' )) + ggtitle('Cadwell')
  }
  p4 = ephysExprPlot(foldy_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T, jitter = T, patchseq = T)
  if (!is.null(p4)){
    p4 = p4 + xlab(paste(gene_name , '(TPM + 1)' )) + ggtitle('Földy')
  }
  p5 = ephysExprPlot(fuzik_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T, jitter = T, patchseq = T)
  if (!is.null(p5)){
    p5 = p5 + xlab(paste(gene_name , '(CPM + 1)' )) + ggtitle('Fuzik')
  }
  
  grid = plot_grid(p2, p3, p4, p5, nrow = 2)
  return(grid)
}

patchSeqCompPlot = function(gene_name, ephys_name, lit_df = discovery_data, aibs_df = jdf, 
                            foldy_df = foldy_joined, cadwell_df = cadwell_joined, fuzik_df = fuzik_joined){
  # lit_df = j
  # aibs_df = jdf
  # foldy_df = foldy_j
  # cadwell_df = cadwell_joined
  
  if (!(gene_name %in% colnames(lit_df))){
    disc_gene_name = 'Hcn3'
    p1 = plot.new()
  } else{
    disc_gene_name = gene_name
    p1 = ephysExprPlot(lit_df, gene_name, ephys_name, add_trend_line=T) + ggtitle('NeuExp/NeuElec')
  }
  
  
  # p2 = ephysExprPlot(aibs_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T)

  
  p2 = ggplot(data = aibs_df, aes_string(x = gene_name, y = ephys_name, size = 'sample_size_weights', weight = 'sample_size_weights', color = 'colors', group = 'broad_type')) + 
    geom_smooth(method = 'lm', aes(group=broad_type, color = broad_type_colors), se = FALSE, alpha = .5, linetype = 2) + 
    geom_smooth(method = 'lm', aes(group=1), se = FALSE, color="grey", alpha = .5, linetype = 1) + 
    geom_point() + 
    scale_color_identity() + 
    geom_point(shape = 1, na.rm=TRUE, aes(size = sample_size_weights, color = broad_type_colors), stroke = 1) +
    scale_x_continuous(trans = 'log10') + annotation_logticks(sides = "b")  + 
    theme(legend.position="none")
  p2 = p2 + xlab(paste(gene_name , '(TPM + 1)' )) + ylab(getEphysNiceName(ephys_name)) + ggtitle('AIBS')
  
  
  p3 = ephysExprPlot(cadwell_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T, jitter = T, patchseq = T)  
  if (!is.null(p3)){
    p3 = p3 + xlab(paste(gene_name , '(TPM + 1)' )) + ggtitle('Cadwell')
  }
  p4 = ephysExprPlot(foldy_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T, jitter = T, patchseq = T)
  if (!is.null(p4)){
    p4 = p4 + xlab(paste(gene_name , '(TPM + 1)' )) + ggtitle('Földy')
  }
  p5 = ephysExprPlot(fuzik_df, gene_name, ephys_name, aibs_flag = T, log_trans_x = T, add_trend_line=T, jitter = T, patchseq = T)
  if (!is.null(p5)){
    p5 = p5 + xlab(paste(gene_name , '(mols + 1)' )) + ggtitle('Fuzik')
  }
  
  grid = plot_grid(p1, p2, plot.new(), p3, p4, p5, nrow = 2)
  return(grid)
}

litAIBSPlot = function(lit_df, aibs_df, gene_name, ephys_name, 
                       show_labels = F, log_trans = T, mult_plot_flag = F, extra_df = F, add_trend_line = T,
                       plot_corr = F, comb_corr_df = NA, remove_y_lab = T) {
  if (plot_corr == F){
    p1 = ephysExprPlot(lit_df, gene_name, ephys_name, show_labels, log_trans, aibs_flag = F, mult_plot_flag, remove_y_lab, add_trend_line = add_trend_line)
    p2 = ephysExprPlot(aibs_df, gene_name, ephys_name, show_labels, log_trans, aibs_flag = T, mult_plot_flag, remove_y_lab, add_trend_line = add_trend_line)
  } else{
    ephys = ephys_name
    disc_corr_row = comb_corr_df %>% dplyr::filter(gene == gene_name, ephys_name == ephys)
    disc_corr_val = round(disc_corr_row$corr, 2)
    disc_corr_pval = signif(10^(-disc_corr_row$fdr), 2)
    aibs_corr_val = round(disc_corr_row$aibs_corrs, 2)

    p1 = ephysExprPlot(lit_df, gene_name, ephys_name, show_labels, log_trans, aibs_flag = F, mult_plot_flag, 
                       remove_y_lab, add_trend_line = add_trend_line, 
                       plot_corr = T, corr_val = disc_corr_val, corr_pval = disc_corr_pval
                       )
    p2 = ephysExprPlot(aibs_df, gene_name, ephys_name, show_labels, log_trans, 
                       aibs_flag = T, mult_plot_flag, remove_y_lab, add_trend_line = add_trend_line,
                       plot_corr = T, corr_val = aibs_corr_val)
  }
  
  # needed for yoking p1 and p2's axes - doesn't look that good
  # p1_ylim = ggplot_build(p1)$layout$panel_ranges[[1]]$y.range
  # p2 = p2 + ylim(p1_ylim)
  
  if (mult_plot_flag){
    POINTSIZE = 2
    TITLESIZE = 10
    TEXTSIZE = 8
    LABELSIZE = 3
  } else{
    POINTSIZE = 4
    TITLESIZE = 16
    TEXTSIZE = 12
    LABELSIZE = 4
  }
  if (is.data.frame(extra_df)){
    limits <- aes(ymax = resp + se, ymin=resp - se)
    p3 = ggplot(extra_df, aes(y=resp, x=trt))
    p3 = p3 + geom_point(stat="identity") + geom_errorbar(limits, width=0.25)
    if (remove_y_lab){
      p3 = p3 + theme(
        axis.title.y = element_blank())
    }
    p3 = p3 + xlab(extra_df$ref[1])
    p3 = p3 + theme(axis.title = element_text(size = TITLESIZE), axis.text = element_text(size=TEXTSIZE))
    
    grid = plot_grid(p1, p2, p3, ncol = 3)
  } else {
    grid = plot_grid(p1, p2, ncol = 2)
  }
}

litAIBSMultGenePlot = function(lit_df, aibs_df, gene_list, ephys_list, show_labels = F, log_trans = T, lit_val_df_list = F, 
                               add_trend_line = F, plot_corr = F, comb_corr_df = NA, grid_columns = 1, remove_y_lab = T) {
  
  plot_list = vector("list", length(gene_list))
  for (i in 1:length(gene_list)){
    if (!is.list(lit_val_df_list)){
      p = litAIBSPlot(lit_df, aibs_df, gene_list[i], ephys_list[i], show_labels, log_trans, mult_plot_flag = T, 
                      add_trend_line = add_trend_line, plot_corr = plot_corr, comb_corr_df = comb_corr_df, remove_y_lab = remove_y_lab)
    } else{
      p = litAIBSPlot(lit_df, aibs_df, gene_list[i], ephys_list[i], show_labels, log_trans, mult_plot_flag = T, 
                      extra_df = lit_val_df_list[[i]], add_trend_line = add_trend_line, plot_corr, comb_corr_df, remove_y_lab = remove_y_lab)
    }
    
    plot_list[[i]] = p
  }
  plot_out = do.call(plot_grid, c(plot_list, ncol = grid_columns, align = 'v'))
}


sampleSummaryPlot = function(df, col_name, sample_names = NULL, sort_names = NULL, 
                             rnaseq_expr_flag = F, take_y_var_log = T) {
  log_vars = c("rin", "tau", "rheo", "maxfreq", "cap")
  
  if (!is.null(sample_names)){
    newdf = df %>% filter(short.name %in% sample_names)
  } else {
    newdf = df
  }
  
  if (!is.null(sort_names)){
    newdf$short.name <- factor(newdf$short.name, levels=sort_names)
    
  }
  
  if (rnaseq_expr_flag == F){
    p = ggplot(newdf, aes_string(x = "short.name", y = col_name, color = "colors"))
  } else{
    p = ggplot(newdf, aes_string(x = "short.name", y = col_name, color = "colors"))
  }
  p = p + geom_violin(na.rm = T, scale = "width", color = "black") 
  #p = p + geom_boxplot(colour = "black", outlier.shape=NA, width = 1) 


  p = p + geom_jitter( na.rm=TRUE, size=2, 
                       position = position_jitter(width = .33), 
                       alpha = .75)
  p = p + scale_colour_identity() 
  #p = p + geom_point(shape = 1,size = 2,colour = "black", na.rm=TRUE, position = position_jitter(width = .5))
  
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
                axis.title.x=element_blank(),
                axis.text.y = element_text(size = 12),
                axis.title.y=element_text(size = 14)
  )
  if ((col_name %in% log_vars) & (take_y_var_log)){
    p = p + scale_y_log10() + annotation_logticks(sides = "l")
  }
  if (rnaseq_expr_flag){
    p = p + scale_y_continuous(trans = "log2") + annotation_logticks(sides = "l")
  }
  p 
}


getGenesFromGO = function(go_term, finalExp){
  tf_probes = list()
  for (i in 1:nrow(finalExp)){
    curr_probe = as.character(finalExp[i, "Probe"])
    GOTerms_w_parents = as.character(GPL1261.an[GPL1261.an$ProbeName == curr_probe,'GOTerms'])
    if (grepl(go_term,GOTerms_w_parents))
      tf_probes = c(tf_probes, curr_probe)
  }
  
  tf_probes = tf_probes %>% unlist
  
  go_group_gene_names = finalExp %>% filter(Probe %in% tf_probes) %>% select(Gene.Symbol) %>% unlist %>% make.names %>% as.vector
  
}

getSigGenesFromGeneList = function(genelist, sig_genes_df_all, sig_genes_df, 
                                   ephys_corr_mat, expr_mat, 
                                   aibs.cormat.spearman, comb_corr_df, min_fdr = .05,
                                   sort_genes_by_sim = T, aibs_corr_match = 0
) {
  
  
  # find genes that are consistent with at least 1 ephys in aibs
  aibs_cons_genes = comb_corr_df %>% filter(cons == 1, abs(aibs_corrs) > aibs_corr_match) %>% select(gene) %>% unlist
  
  # try to aggregate all significant genes across properties
  ephys_gene_mat = sig_genes_df_all %>% 
    select(gene, ephys_name, corr) %>% 
    filter(gene %in% genelist) %>% 
    #filter(gene %in% sig_genes_df$gene) %>% # does filtering for at least 1 gene in sig list
    filter(gene %in% aibs_cons_genes, gene %in% sig_genes_df$gene) %>%
    spread(ephys_name, corr)
  rownames(ephys_gene_mat) = ephys_gene_mat$gene
  ephys_gene_mat$gene = NULL
  
  
  m = ephys_gene_mat %>% as.matrix
  
  
  temp_ephys_names = 
    rownames(ephys_corr_mat)
  
  # source functions for sorting corr mat
  source('~/ephys_analysis/sortCorrMat.R')
  
  ## plot correlations among ephys data
  sorted_ephys_col_names = colnames(ephys_corr_mat[,rownames(ephys_corr_mat)  %in%  colnames(m)])
  sorted_ephys_col_names = sorted_ephys_col_names [sorted_ephys_col_names %in% colnames(m)]
  
  # calculate pairwise expr correlation matrix for sorting purposes
  if (sort_genes_by_sim){
    sorted_gene_names = rownames(reorder_cormat(cor(expr_mat[,rownames(ephys_gene_mat )], , use="pairwise.complete.obs")))
  } else{
    sorted_gene_names = rev(rownames(ephys_gene_mat ))
  }
  
  # Melt the correlation matrix
  melted_ephys_gene_mat <- melt(t(m[sorted_gene_names, sorted_ephys_col_names]), na.rm = TRUE)
  if (sort_genes_by_sim == F){
    melted_ephys_gene_mat$Var2 = factor(melted_ephys_gene_mat$Var2, levels = rownames(ephys_gene_mat ))
  }
  
  melted_ephys_gene_mat$fdr = ''
  for (i in 1:nrow(melted_ephys_gene_mat)){
    ephys = melted_ephys_gene_mat[i,'Var1']
    gene_name = melted_ephys_gene_mat[i,'Var2']
    fdr_val = sig_genes_df_all %>% filter(gene == gene_name, ephys_name == ephys) %>% select(fdr)
    
    
    aibs_text = ''
    cons = 0
    nt_adj_pval = 1
    
    # figure out aibs code
    if (fdr_val > -log10(.1)){
      if (gene_name %in% colnames(aibs.cormat.spearman) & ephys %in% rownames(aibs.cormat.spearman)) {
        aibs_row = comb_corr_df %>% filter(gene == gene_name, ephys_name == ephys)
        if (nrow(aibs_row) > 0){
          aibs_corr_val = aibs_row$aibs_corrs
          brainmash_corr_val = aibs_row$corr
          nt_adj_pval = sig_genes_df %>% filter(gene == gene_name, ephys_name == ephys_name) %>% select(nt_adj_pval)
          #print (nt_adj_pval)
          if ( (abs(aibs_corr_val) > aibs_corr_match) & (aibs_row$cons == 1)){
            #if ( (abs(brainmash_corr_val) - abs(aibs_corr_val)) < .3){
            cons = 1
          }
          #cons = aibs_row$cons
          if (cons == 1){
            aibs_text = ''
          } else{
            aibs_text = '/'
          }
        }
      } else{
        aibs_text = '?'
      }
#       if (nt_adj_pval > .1){
#         aibs_text = paste0('x',aibs_text)
#       }
    }

    if (fdr_val > -log10(.1) & fdr_val < -log10(.05)){
      melted_ephys_gene_mat[i,'fdr'] = paste0('·',aibs_text)
    } else if (fdr_val > -log10(.05) & fdr_val < -log10(.01)){
      melted_ephys_gene_mat[i,'fdr'] = paste0('*',aibs_text)
    } else if (fdr_val > -log10(.01)) {
      melted_ephys_gene_mat[i,'fdr'] = paste0('**',aibs_text)
    }
    
    
  }
  
  response_vars_plot_names = c('Rin', 'Vrest', 'APthr', 'APamp', 'APhw', 'Tau', 'AHPamp', 'Rheo', 'FRmax', 'Cm', 'SFA')

  melted_ephys_gene_mat$Var1 = plyr::mapvalues(melted_ephys_gene_mat$Var1, from = response_vars, to = response_vars_plot_names)

  
  # Create a ggheatmap
  ephys_gene_sig_heatmap <- ggplot(melted_ephys_gene_mat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorr.") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0, 
                                     size = 8, hjust = 1),
          axis.text.y = element_text(size = 8, hjust = 1))+
    coord_fixed()+ 
    scale_x_discrete(position = "top") + 
    geom_text(aes(Var2, Var1, label = fdr), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x=element_text(vjust=1, hjust=0, color = "black"),
      axis.text.y=element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "vertical", legend.position="right",
      legend.title =element_text(size=10),
      legend.text = element_text(size=8)
      ) 

  # Print the heatmap
  ephys_gene_sig_heatmap
}

ageDepthPlot = function(exprname, ephysname){
  
  ephys_name = ephysname
  gene_name = exprname
  
  log_vars = c("rin", "tau", "rheo", "maxfreq", "cap")
    

  
  p1 = ggplot(pv_data_expr, aes_string(x = 'age', y = gene_name, color = "colors"))
  p1 = p1 + geom_point(na.rm=TRUE, size=2) + geom_smooth(se =F)
  p1 = p1 +  scale_colour_identity() + xlim(5, 40)

  p2 = ggplot(pv_data_ephys, aes_string(x = 'AnimalAge', y = ephys_name, color = "colors"))
  p2 = p2 + geom_point(na.rm=TRUE, size=2)  + geom_smooth(se =F)
  p2 = p2 +  scale_colour_identity()  +  xlim(5, 40)
  #p2 = p2 + ylim(-80, -50)
  
  
  p3 = ggplot(purk_data_expr, aes_string(x = 'age', y = gene_name, color = "colors"))
  p3 = p3 + geom_point(na.rm=TRUE, size=2)  + geom_smooth(se =F)
  p3 = p3 +  scale_colour_identity()  + xlim(2, 60)
  
  p4 = ggplot(purk_data_ephys, aes_string(x = 'AnimalAge', y = ephys_name, color = "colors"))
  p4 = p4 + geom_point(na.rm=TRUE, size=2)  + geom_smooth(se =F)
  p4 = p4 +  scale_colour_identity()  
  
  if ((ephys_name %in% log_vars) ){
    p2 = p2 + scale_y_log10() + annotation_logticks(sides = "l")
    p4 = p4 + scale_y_log10() + annotation_logticks(sides = "l")
  }
  #p4 = p4 + ylim(-80, -50)
  plot_grid(p1, p3, p2, p4, ncol=2, align = 'v')

}
