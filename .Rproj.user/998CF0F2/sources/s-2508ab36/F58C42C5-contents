source('R/calcCellContam.R') # used for getValidGenes function
source('R/getGoodMarkers.R')
source('data-raw/markers/NeuroExpressoMarkers.R')

### define ON cell type markers - markers that are expressed in the cell subtype of interest

print("Defining ON cell type markers based on analysis of Tasic and Zeisel datasets")
on_cell_types = c('Ndnf', 'Sncg', 'Pvalb', 'Pyramidal')

# for inhibitory and pyramidal cells, what other cell types are we comparing these to?
off_broad_types_inh = c('Astrocyte', 'Microglia', 'Pyramidal', 'Oligodendrocyte', 'OPC', 'Endothelial')
off_broad_types_pyr = c('Astrocyte', 'Microglia', 'Inhibitory', 'Oligodendrocyte', 'OPC', 'Endothelial')

on_cell_compare_types = list(Ndnf = off_broad_types_inh, Sncg = off_broad_types_inh, Pvalb = off_broad_types_inh,
                             Pyramidal = off_broad_types_pyr)
examplar_patch_seq_datasets = list(Ndnf = 'cadwell', Sncg = 'foldy', Pvalb = 'foldy', Pyramidal = 'fuzik')

# workflow for on markers
# get initial set of markers based on Tasic
# filter markers for goodness based on expression in Tasic and Zeisel
# sort markers based expression levels in relevant patch-seq dataset

RATIO_THRESHOLD = 10
MEAN_ON_CELL_EXPR_THRESHOLD = 100

ON_MARKER_THRESHOLD_TASIC = 10 + 1 #TPM +1
ON_MARKER_THRESHOLD_ZEISEL = .01 # UMI counts
MIN_CELLS_EXPRESSING_RATIO_TASIC = .75
MIN_CELLS_EXPRESSING_RATIO_ZEISEL = .5

onCellMarkerGenes = lapply(on_cell_types, function(cell_type_name){
  curr_on_cell_type = cell_type_name


  # ndnf_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% ndnf_cells, aibs_gene_names] )
  # pvalb_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% pvalb_cells, aibs_gene_names] )

  on_cell_inds = which(aibsExprDataDf$norm_sub_types == curr_on_cell_type)
  on_mean_expr = aibsExprDataDf[on_cell_inds, aibs_gene_names] %>% colMeans() # note that aibs data is in TPM + 1

  off_cell_inds = which(aibsExprDataDf$norm_broad_types %in% on_cell_compare_types[[curr_on_cell_type]])
  off_mean_expr = aibsExprDataDf[off_cell_inds, aibs_gene_names] %>% colMeans() # note that aibs data is in TPM + 1

  marker_ratio = on_mean_expr / off_mean_expr
  # defines initial list of markers based on tasic
  initial_markers = marker_ratio[marker_ratio > RATIO_THRESHOLD & on_mean_expr > MEAN_ON_CELL_EXPR_THRESHOLD] %>% names %>% unlist %>% as.character()

  filtered_markers_tasic = getGoodMarkersFxn(curr_on_cell_type, initial_markers, 'tasic',
                                          MARKER_THRESHOLD = ON_MARKER_THRESHOLD_TASIC, MIN_CELLS_EXPRESSING_RATIO = MIN_CELLS_EXPRESSING_RATIO_TASIC)
  filtered_markers_zeisel = getGoodMarkersFxn(curr_on_cell_type, getValidGenes(initial_markers, zeisel_gene_names), 'zeisel',
                                           MARKER_THRESHOLD = ON_MARKER_THRESHOLD_ZEISEL, MIN_CELLS_EXPRESSING_RATIO =MIN_CELLS_EXPRESSING_RATIO_ZEISEL)
  consensus_on_genes = intersect(filtered_markers_tasic, filtered_markers_zeisel)

  # now sort genes by expression levels in example patch seq datasets

  patch_seq_expr = patch_seq_datasets[[examplar_patch_seq_datasets[[curr_on_cell_type]]]]$joined_df
  patch_seq_rel_cells = patch_seq_datasets[[examplar_patch_seq_datasets[[curr_on_cell_type]]]]$joined_df$contam_type == curr_on_cell_type
  patch_seq_cell_expr = patch_seq_expr[patch_seq_rel_cells, consensus_on_genes]
  sorted_marker_expr =  patch_seq_cell_expr %>% colMeans %>% sort(., decreasing = T) %>% names() %>% unlist %>% as.character

  return(sorted_marker_expr)
})
names(onCellMarkerGenes) = on_cell_types


### now define list of off cell type markers
print("Defining OFF cell type markers based on analysis of Tasic and Zeisel datasets")

# get markers from neuroexpresso
fullMarkerList = neuroExpressoMarkers$singleCell

# markers from broad cell types
useCellTypeList = c('Astrocyte', 'Microglia', 'Oligodendrocyte', 'Pyramidal', 'Endothelial', 'Oligodendrocyte precursors')
fullMarkerList = fullMarkerList[names(fullMarkerList) %in% useCellTypeList]
# names(fullMarkerList)[which(names(fullMarkerList) == "FS Basket (G42)")] = 'Pvalb'
names(fullMarkerList)[which(names(fullMarkerList) == 'Oligodendrocyte precursors')] = 'OPC'

# define some initial pan-gabargic markers
gaba_cells = aibsExprDataDf %>% select(one_of('broad_type', 'sample_title')) %>% filter(str_detect(broad_type, 'GABA')) %>% select(sample_title) %>% unlist
non_gaba_cells = aibsExprDataDf %>% select(one_of('broad_type', 'sample_title')) %>% filter(!str_detect(broad_type, 'GABA')) %>% select(sample_title) %>% unlist
gaba_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% gaba_cells, aibs_gene_names] )

non_inh_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% non_gaba_cells, aibs_gene_names] )

mydf = data.frame(non_inh = non_inh_expr_profile,
                  gaba = gaba_expr_profile) %>% tibble::rownames_to_column(var = "gene")

gaba_markers = mydf %>% mutate(ratio = gaba / non_inh) %>%
  filter(ratio > 10, gaba > 100) %>% arrange(-gaba)  %>% select(gene) %>% unlist() %>% as.character
# gaba_markers = getGoodMarkers(gaba_cells, gaba_markers)

fullMarkerList$Inhibitory = gaba_markers


off_cell_types = c('Astrocyte', 'Microglia', 'Inhibitory', 'Oligodendrocyte', 'OPC', 'Endothelial', 'Pyramidal')

offCellMarkerGenes = lapply(off_cell_types, function(off_cell_type){
  curr_off_cell_type = off_cell_type

  # step 1: are markers expressed at reasonable levels in cells they should be expressed in?

  tasic_off_markers = getGoodMarkersFxn(curr_off_cell_type, getValidGenes(fullMarkerList[[curr_off_cell_type]], aibs_gene_names), 'tasic',
                                     MARKER_THRESHOLD = ON_MARKER_THRESHOLD_TASIC, MIN_CELLS_EXPRESSING_RATIO = .5, use_subtype = F)
  zeisel_off_markers = getGoodMarkersFxn(curr_off_cell_type, getValidGenes(fullMarkerList[[curr_off_cell_type]], zeisel_gene_names), 'zeisel',
                                      MARKER_THRESHOLD = ON_MARKER_THRESHOLD_ZEISEL, MIN_CELLS_EXPRESSING_RATIO = .5, use_subtype = F)


  ok_off_markers = intersect(tasic_off_markers, zeisel_off_markers)
  # print(ok_off_markers)
  # return(ok_off_markers)
  # step 2: now go through every on cell type and make sure markers aren't really expressed in that on cell type

  bad_markers = lapply(on_cell_types, function(on_cell_type){
    # don't compare pyramidal cell markers to pyramidal cells
    if(on_cell_type == curr_off_cell_type){
      return(NULL)
    }
    # don't compare inhibitory cell markers to inhibitory cells
    if((curr_off_cell_type == 'Inhibitory') & (on_cell_type %in% c('Pvalb', 'Sncg', 'Ndnf'))){
      return(NULL)
    }

    bad_markers_tasic = getGoodMarkersFxn(on_cell_type, ok_off_markers, 'tasic',
                                       MARKER_THRESHOLD = 10 + 1, MIN_CELLS_EXPRESSING_RATIO = .33, marker_direction = 'positive', use_subtype = T)

    bad_markers_zeisel = getGoodMarkersFxn(on_cell_type, ok_off_markers, 'zeisel',
                                       MARKER_THRESHOLD = 2, MIN_CELLS_EXPRESSING_RATIO = .33, marker_direction = 'positive', use_subtype = T)

    curr_bad_markers = union(bad_markers_tasic, bad_markers_zeisel)
    # curr_bad_markers = bad_markers_tasic
    # print(on_cell_type)
    return(curr_bad_markers)

  })
  # print(off_cell_type)
  # print(bad_markers %>% unlist %>% unique())
  names(bad_markers) = on_cell_types
  bad_markers_flat = bad_markers %>% unlist %>% unique()
  ok_off_markers_filtered = setdiff(ok_off_markers, bad_markers_flat)

  return(ok_off_markers_filtered)
})
names(offCellMarkerGenes) = off_cell_types


offCellMarkerGenes = lapply(names(offCellMarkerGenes), function(cell_type){
  patch_seq_expr = patch_seq_datasets[['cadwell']]$joined_df
  patch_seq_rel_cells = patch_seq_datasets[['cadwell']]$joined_df$contam_type == 'Ndnf'
  patch_seq_cell_expr = patch_seq_expr[patch_seq_rel_cells, offCellMarkerGenes[[cell_type]]]
  sorted_marker_expr =  patch_seq_cell_expr %>% colMeans %>% sort(., decreasing = T) %>% names() %>% unlist %>% as.character
  return(sorted_marker_expr)
})
names(offCellMarkerGenes) = off_cell_types

fullMarkerList[names(offCellMarkerGenes)] = offCellMarkerGenes

names(onCellMarkerGenes) = paste0(names(onCellMarkerGenes), '_on')

fullMarkerList[names(onCellMarkerGenes)] = onCellMarkerGenes

### output list of markers to a file

marker_df = do.call(rbind, lapply(seq_along(fullMarkerList), function(i){
  data.frame(CLUSTER=names(fullMarkerList)[i], fullMarkerList[[i]])
}))

colnames(marker_df) = c('cell_type', 'gene')

write.csv(marker_df, 'analysis/tables/markerGeneList.csv')

markers= fullMarkerList

devtools::use_data(markers, overwrite = T)


