### load lists of markers from neuroexpresso and other tools

# get markers from neuroexpresso

source('data-raw/markers/NeuroExpressoMarkers.R')
# produces a variable named neuroExpressoMarkers

### calculate marker lists for specific cell types based on Tasic 2016 data

print('Calculating extra markers from Tasic for Ndnf and Sncg cells')

print('Calculating average gene expression for Tasic 2016 subtypes -- this takes a while')
aibs_primary = aibsExprDataDf %>% group_by(primary_type) %>% summarize_at(aibs_gene_names, mean)
rownames(aibs_primary) = aibs_primary$primary_type

aibs_primary = left_join(aibs_primary, aibsExprDataDf %>% select(primary_type, broad_type) %>% distinct(primary_type, broad_type), by = 'primary_type') %>% select(primary_type, broad_type, everything())
rownames(aibs_primary) = aibs_primary$primary_type


aibs_primary_trans = aibs_primary[, aibs_gene_names] %>% t() %>% as.data.frame %>% tibble::rownames_to_column(var = 'gene')
#rownames(exprdf_trans) = exprdf_trans$primary_type
colnames(aibs_primary_trans)= c('gene', aibs_primary$primary_type %>% unlist %>% as.character() %>% make.names())
rownames(aibs_primary_trans) = aibs_primary_trans$gene
#exprdf_trans = exprdf_trans[, 2:ncol(exprdf_trans)]

aibs_exc_names = aibs_primary %>% select(one_of(c('primary_type', 'broad_type'))) %>% filter(str_detect(broad_type, 'Gluta')) %>% select(primary_type) %>% unlist %>% make.names()

ndnf_cells = aibsExprDataDf %>% select(broad_type, primary_type, sample_title) %>% filter(str_detect(primary_type, 'Ndnf')) %>% select(sample_title) %>% unlist
sncg_cells = aibsExprDataDf %>% select(broad_type, primary_type, sample_title) %>% filter(str_detect(primary_type, 'Sncg')) %>% select(sample_title) %>% unlist
pvalb_cells = aibsExprDataDf %>% select(broad_type, primary_type, sample_title) %>% filter(str_detect(primary_type, 'Pvalb')) %>% select(sample_title) %>% unlist


sncg_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% sncg_cells, aibs_gene_names] )

ndnf_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% ndnf_cells, aibs_gene_names] )
pvalb_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% pvalb_cells, aibs_gene_names] )


exc_cells = aibsExprDataDf %>% select(broad_type, primary_type, sample_title) %>% filter(str_detect(broad_type, 'Gluta')) %>% select(sample_title) %>% unlist
gaba_cells = aibs_meta %>% filter(str_detect(broad_type, 'GABA')) %>% select(sample_title) %>% unlist
non_gaba_cells = aibs_meta %>% filter(!str_detect(broad_type, 'GABA')) %>% select(sample_title) %>% unlist
non_exc_cells = aibs_meta %>% filter(!str_detect(broad_type, 'Gluta')) %>% select(sample_title) %>% unlist

aibs_meta = aibsExprDataDf %>% select(broad_type, primary_type, sample_title) 

exc_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% exc_cells, aibs_gene_names] )
gaba_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% gaba_cells, aibs_gene_names] )

non_inh_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% non_gaba_cells, aibs_gene_names] )
non_exc_expr_profile = colMeans(aibsExprDataDf[aibsExprDataDf$sample_title %in% non_exc_cells, aibs_gene_names] )


mydf = data.frame(ndnf = ndnf_expr_profile, exc = exc_expr_profile, non_inh = non_inh_expr_profile, sncg = sncg_expr_profile,
                  non_exc = non_exc_expr_profile, 
                  gaba = gaba_expr_profile,
                  pvalb = pvalb_expr_profile) %>% tibble::rownames_to_column(var = "gene")



sortMarkers = function(aibs_cell_names, markers){
  sorted_markers = aibsExprDataDf[aibsExprDataDf$sample_title %in% aibs_cell_names, markers] %>% 
    summarize_all(median) %>% t() %>% 
    as.data.frame() %>% tibble::rownames_to_column(var = 'gene') %>%
    rename(median_expr = V1) %>% arrange(median_expr) %>% 
    select(gene) %>% unlist %>% as.character
  return(sorted_markers)
}

ndnf_markers = mydf %>% mutate(ratio = ndnf / non_inh) %>% 
  filter(ratio > 10, ndnf > 100) %>% arrange(-ndnf)  %>% select(gene) %>% unlist() %>% as.character
ndnf_markers = getGoodMarkers(ndnf_cells, ndnf_markers)

sncg_markers = mydf %>% mutate(ratio = sncg / non_inh) %>% 
  filter(ratio > 10, sncg > 100) %>% arrange(-sncg)  %>% select(gene) %>% unlist() %>% as.character
sncg_markers = getGoodMarkers(sncg_cells, sncg_markers)
pvalb_markers = mydf %>% mutate(ratio = pvalb / non_inh) %>% 
  filter(ratio > 10, pvalb > 100) %>% arrange(-pvalb)  %>% select(gene) %>% unlist() %>% as.character
pvalb_markers = getGoodMarkers(pvalb_cells, pvalb_markers)

pyr_markers = mydf %>% mutate(ratio = exc / non_exc) %>% 
  filter(ratio > 10, exc > 100) %>% arrange(-exc)  %>% select(gene) %>% unlist() %>% as.character
pyr_markers = getGoodMarkers(exc_cells, pyr_markers)
gaba_markers = mydf %>% mutate(ratio = gaba / non_inh) %>% 
  filter(ratio > 10, gaba > 100) %>% arrange(-gaba)  %>% select(gene) %>% unlist() %>% as.character
gaba_markers = getGoodMarkers(gaba_cells, gaba_markers)


ndnf_markers = aibs_primary_trans %>% rowwise() %>% mutate_at(ndnf = mean(c(`Ndnf.Car4`, `Ndnf.Cxcl14`)), 
                                             l23 = mean(paste),
                                             ratio = ndnf / l23) %>% 
  filter(ratio > 50, ndnf > 200) %>% arrange(-ndnf) %>% select(gene) %>% unlist() %>% as.character

ndnf_markers = aibs_primary_trans %>% mutate(ndnf = rowMeans(cbind(`Ndnf.Car4`, `Ndnf.Cxcl14`)), 
                                             l23 = rowMeans(aibs_exc_names),
                                             ratio = ndnf / l23) %>% 
  filter(ratio > 50, ndnf > 200) %>% arrange(-ndnf) %>% select(gene) %>% unlist() %>% as.character

ndnf_markers = aibs_primary_trans %>% mutate(ndnf = rowMeans(cbind(`Ndnf.Car4`, `Ndnf.Cxcl14`)), 
                                       l23 = rowMeans(cbind(`L2.Ngb`,`L2.3.Ptgs2`)),
                                       ratio = ndnf / l23) %>% 
  filter(ratio > 50, ndnf > 200) %>% arrange(-ndnf) %>% select(gene) %>% unlist() %>% as.character

sncg_markers = aibs_primary_trans %>% mutate(sncg = rowMeans(cbind(`Sncg`)), 
                                       l23 = rowMeans(cbind(`L2.Ngb`,`L2.3.Ptgs2`)),
                                       ratio = sncg / l23) %>% 
  filter(ratio > 50, sncg > 200) %>% arrange(-sncg) %>% select(gene) %>% unlist() %>% as.character

# filter out sparsely 
nnz_gene_counts = aibsExprDataDf %>% select(one_of(c(sncg_markers, 'primary_type'))) %>% 
  filter(primary_type == 'Sncg') %>% select(-primary_type) %>% summarize_all(funs(sum(. >= 10)))
sncg_markers = sncg_markers[which(nnz_gene_counts >=8)]


fullMarkerList = neuroExpressoMarkers$singleCell
# fullMarkerList = neuroExpressoMarkers$microarray$Cortex

# filter out gabaergic neuroexpresso types
useCellTypeList = c('Astrocyte', 'Microglia', 'Oligodendrocyte', 'Pyramidal', 'FS Basket (G42)', 'Endothelial', 'Oligodendrocyte precursors')
fullMarkerList = fullMarkerList[names(fullMarkerList) %in% useCellTypeList]
names(fullMarkerList)[which(names(fullMarkerList) == "FS Basket (G42)")] = 'Pvalb'
names(fullMarkerList)[which(names(fullMarkerList) == 'Oligodendrocyte precursors')] = 'OPC'
fullMarkerList$Inhibitory = gaba_markers

onCellMarkers = list(Ndnf = ndnf_markers, Sncg = sncg_markers, Pvalb = pvalb_markers, Pyramidal = pyr_markers)


fullMarkerList$Ndnf = ndnf_markers
fullMarkerList$Sncg = sncg_markers
fullMarkerList$Pvalb = pvalb_markers
# fullMarkerList$Pyramidal = pyr_markers



# try filtering tasic cell types
TPM_THRESHOLD = 5
compare_marker_types =  c('Astrocyte', 'Microglia', 'Oligodendrocyte', 'Endothelial', 'OPC')

offCellMarkers = lapply(compare_marker_types, function(marker_type){
  curr_marker_names = intersect(fullMarkerList[marker_type] %>% unlist %>% as.character(), aibs_gene_names)
  gaba_cell_marker_expr = gaba_expr_profile[curr_marker_names]
  gaba_cell_valid_markers = gaba_cell_marker_expr[gaba_cell_marker_expr < TPM_THRESHOLD + 1] %>% names
  
  
  exc_cell_marker_expr = exc_expr_profile[curr_marker_names]
  exc_cell_valid_markers = exc_cell_marker_expr[exc_cell_marker_expr < TPM_THRESHOLD + 1] %>% names
  
  updated_markers = intersect(gaba_cell_valid_markers, exc_cell_valid_markers)
  return(updated_markers)
})
names(offCellMarkers) = compare_marker_types

# filter pyramidal markers
off_pyramidal_markers = fullMarkerList$Pyramidal

# filter out pyramidal markers expressed too highly in inhibitory cells
gaba_cell_marker_expr = gaba_expr_profile[getValidGenes(off_pyramidal_markers, aibs_gene_names)]
gaba_cell_valid_markers = gaba_cell_marker_expr[gaba_cell_marker_expr < TPM_THRESHOLD + 1] %>% names


offCellMarkers$Pyramidal = sort(-exc_expr_profile[gaba_cell_valid_markers]) %>% names

offCellMarkers$Pyramidal =  colSums(cadwell_expr[offCellMarkers$Pyramidal]) %>% sort(,decreasing = T) %>% names


fullMarkerList[names(onCellMarkers)] = onCellMarkers
fullMarkerList[names(offCellMarkers)] = offCellMarkers

# fullMarkerList$Sncg = neuroExpressoMarkers$singleCell$`VIPReln (G30)`

save(aibs_primary, ndnf_markers, sncg_markers, file = 'data/tasic_extra_markers.rda')
save(fullMarkerList, file = 'data/fullMarkerList.rda') # not filtered yet for zeisel expression



