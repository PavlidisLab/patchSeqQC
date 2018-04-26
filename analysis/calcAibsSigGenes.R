source('R/corr_funcs.R')

print('Loading data frame with aibs / Tasic ephys and expression data summarized to the cre-line and layer level')
aibs_cre_joined = readRDS('data/aibs_cre_expr_ephys.rda')


## calculate correlations between ephys and expression
print('Calculating pairwise correlations between gene expr and ephys features for aibs / Tasic data')

response_vars = c("rin", "rmp","apthr","apamp",
                  "aphw", "tau", "ahpamp", "rheo","maxfreq", "cap", "adratio")
log_vars = c("rin", "tau", "aphw", "rheo", "maxfreq", "cap")

inh_cell_types = which(aibs_cre_joined$broad_type == 'Inh')

use_cell_types = which(aibs_cre_joined$broad_type %in% c('Exc', 'Inh'))
aibs_expr_mat = aibs_cre_joined[use_cell_types,aibs_gene_names] %>% as.matrix
aibs_ephys_mat = aibs_cre_joined[use_cell_types,response_vars] %>% as.matrix

aibs_expr_mat_sd = apply(aibs_expr_mat, 2, sd)

aibs_ephys_corrs = calculateExprEphysWeightedCorr(aibs_expr_mat,  aibs_ephys_mat, weights = aibs_cre_joined[use_cell_types, "sample_size_weights"],
                                                  ephys_var_names = response_vars, log_ephys_vars = log_vars)
