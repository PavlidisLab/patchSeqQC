### load dissociated cell reference data from .Rda frames stored on disk
### contact stripathy for csv files if required

aibsExprDataDf = readRDS('data/tasic_expression.rda')
aibs_gene_names = readRDS('data/tasic_gene_names.rda')
expr_cell_colors = readRDS('data/tasic_expr_cell_colors.rda')
aibsColors = readRDS('data/aibs_colors.rda')

zeiselExprDataDf = readRDS('data/zeisel_expression.rda')
zeisel_gene_names = readRDS('data/zeisel_gene_names.rda')
# 
# aibs_cre_joined = readRDS('data/aibs_cre_expr_ephys.rda')
# 
# source('data-raw/markers/NeuroExpressoMarkers.R')
# fullMarkerList = readRDS('data/full_marker_list.rda')
# 
# 
# patch_seq_datasets = readRDS('data/patch_seq_datasets.rda')
# load("~/ephys_analysis/patchSeqQC/data/tasic_extra_markers.Rdata")
# aibsEphys = readRDS('data/aibs_ephys.rda')

