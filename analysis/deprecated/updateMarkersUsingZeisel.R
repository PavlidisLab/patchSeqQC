## generate final list of markers that we'll use in plots and analysis

### compare expression of cell markers in zeisel cell types - if weird, then remove the gene from marker list

# Int15 - Ndnf Car4
# Int12 - Ndnf Cxcl14
# Int6 - Sncg

print('Filtering NeuroExpresso single cell markers based on expression in Zeisel')

zeisel_compare_cell_groups = c('Int15', 'Int12', 'Int6') # based on Maggies paper showing that these are the 
compare_marker_types =  c('Astrocyte', 'Microglia', 'Oligodendrocyte', 'Pyramidal', 'Endothelial', 'Oligodendrocyte precursors')

bad_gene_list = c()
zeisel_expr_threshold = .5
o = lapply(zeisel_compare_cell_groups, function(zeisel_cell_type){
  o = lapply(compare_marker_types, function(marker_type){
    marker_genes = fullMarkerList[[marker_type]]
    marker_expr = zeisel_grps_mean[marker_genes, zeisel_cell_type]
    curr_bad_genes = names(marker_expr[marker_expr > zeisel_expr_threshold])
    bad_gene_list = c(bad_gene_list, curr_bad_genes)
    return(bad_gene_list)
  })
  return(o)
})

bad_gene_list = unlist(o) %>% unique()

for (i in 1:length(compare_marker_types)){
  fullMarkerList[[compare_marker_types[i]]] = setdiff(fullMarkerList[[compare_marker_types[i]]], bad_gene_list)
}

### filter Ndnf and Sncg genes taht aren't highly expressed in zeisel

bad_ndnf_genes = zeisel_grps_mean[fullMarkerList$Ndnf, c('Int15', 'Int12')] %>% tibble::rownames_to_column(var = 'gene') %>%  filter(Int15 < zeisel_expr_threshold | Int12 < zeisel_expr_threshold) %>% select(gene) %>% unlist
fullMarkerList$Ndnf = setdiff(fullMarkerList$Ndnf , bad_ndnf_genes)

bad_sncg_genes = zeisel_grps_mean[fullMarkerList$Sncg, c('Int6', 'Int12')] %>% tibble::rownames_to_column(var = 'gene') %>%  filter(Int6 < zeisel_expr_threshold) %>% select(gene) %>% unlist
fullMarkerList$Sncg = setdiff(fullMarkerList$Sncg , bad_sncg_genes)

### now update order of markers by their expression magnitude in zeisel
fullMarkerList$Pyramidal = names(zeisel_grps_mean[fullMarkerList$Pyramidal, 'S1PyrL23'] %>% sort(decreasing = T))


saveRDS(fullMarkerList, file = 'data/full_marker_list.rda')

# 
# ### try updating gene lists to be more specific across broad types
# 
# match_cell_types = list(Astrocyte = c("Astro1", "Astro2"), Ependymal = c("Epend"), 
#                         Pyramidal = c("S1PyrL4",  "ClauPyr",  "S1PyrL5",  "S1PyrL23", "S1PyrDL",  "S1PyrL5a", "S1PyrL6b", "S1PyrL6", 
#                                       "CA1Pyr1", "CA1Pyr2", "CA1PyrInt", "CA2Pyr2", "CA1Pyr1", "SubPyr"),
#                         Oligodendrocyte = c("Oligo1",   "Oligo3", "Oligo4" ,  "Oligo2" ,  "Oligo6" ,  "Oligo5"),
#                         Microglia = c("Mgl1", "Mgl2")
#                         )
# 
# zeisel_expr_threshold = 1
# 
# for (i in 1:length(match_cell_types)){
#   
#   curr_cell_type = match_cell_types[i] %>% names()
#   print(curr_cell_type)
#   zeisel_compare_cell_groups = setdiff(zeisel_cell_types, match_cell_types[[i]])
#   
#   bad_gene_list = c()
# 
#   o = lapply(zeisel_compare_cell_groups, function(zeisel_cell_type){
#     curr_compare_marker_types = compare_marker_types[!compare_marker_types %in% curr_cell_type]
#   
#     marker_genes = fullMarkerList[[curr_cell_type]]
#     marker_expr = zeisel_grps_mean[marker_genes, zeisel_cell_type]
#     curr_bad_genes = names(marker_expr[marker_expr > zeisel_expr_threshold])
#     bad_gene_list = c(bad_gene_list, curr_bad_genes)
#     return(bad_gene_list)
#   })
# 
#   bad_gene_list = unlist(o) %>% unique()
#   
#   # print(curr_cell_type)
#   # print(bad_gene_list)
#   # print(fullMarkerList[[curr_cell_type]])
#   
#   updated_markers = setdiff(fullMarkerList[[curr_cell_type]], bad_gene_list)
#   
#   if (length(match_cell_types[[i]]) > 1)
#     updated_markers = zeisel_grps_mean[updated_markers, match_cell_types[[i]]] %>% rowMeans(na.rm = T)  %>% sort(decreasing = T) %>% names
#   else
#     updated_markers = zeisel_grps_mean[updated_markers, match_cell_types[[i]]] %>% sort(decreasing = T) %>% names
#   
#   print(updated_markers)
#   
#   if (length(updated_markers > 10)){
#     fullMarkerList[[curr_cell_type]] = updated_markers
#   }
# 
# }
# fullMarkerList = fullMarkerList[setdiff(names(fullMarkerList), 'OPC')]
# # fullMarkerList$OPC = c('Neu4', 'Vcan')
# 
# print(setdiff(fullMarkerList$Astrocyte, bad_gene_list))
# 
# zeisel_compare_cell_groups = setdiff(zeisel_cell_types,c("Astro1", "Astro2"))
# 
# # zeisel_compare_cell_groups = c('Mgl1', 'Mgl2', "Oligo1") # based on Maggies paper showing that these are the 
# compare_marker_types =  c('Astrocyte')
# 
# bad_gene_list = c()
# zeisel_expr_threshold = .5
# o = lapply(zeisel_compare_cell_groups, function(zeisel_cell_type){
#   o = lapply(compare_marker_types, function(marker_type){
#     print(marker_type)
#     marker_genes = fullMarkerList[[marker_type]]
#     marker_expr = zeisel_grps_mean[marker_genes, zeisel_cell_type]
#     curr_bad_genes = names(marker_expr[marker_expr > zeisel_expr_threshold])
#     bad_gene_list = c(bad_gene_list, curr_bad_genes)
#     return(bad_gene_list)
#   })
#   return(o)
# })
# 
# bad_gene_list = unlist(o) %>% unique()
# 
# newMarkers = setdiff(fullMarkerList$Astrocyte, bad_gene_list)
# newMarkers = intersect(newMarkers, zeisel_gene_names)




