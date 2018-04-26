### load lists of markers from neuroexpresso and other tools

# get markers from neuroexpresso

source('data-raw/markers/NeuroExpressoMarkers.R')
# produces a variable named neuroExpressoMarkers

### calculate marker lists for specific cell types based on Tasic 2016 data

aibs_primary = aibsExprDataDf %>% group_by(primary_type) %>% summarize_at(aibs_gene_names, mean)
rownames(aibs_primary) = aibs_primary$primary_type

aibs_primary = left_join(aibs_primary, aibsExprDataDf %>% select(primary_type, broad_type) %>% distinct(primary_type, broad_type), by = 'primary_type') %>% select(primary_type, broad_type, everything())
rownames(aibs_primary) = aibs_primary$primary_type


aibs_primary_trans = aibs_primary[, aibs_gene_names] %>% t() %>% as.data.frame %>% tibble::rownames_to_column(var = 'gene')
#rownames(exprdf_trans) = exprdf_trans$primary_type
colnames(aibs_primary_trans)= c('gene', aibs_primary$primary_type %>% unlist %>% as.character() %>% make.names())
rownames(aibs_primary_trans) = aibs_primary_trans$gene
#exprdf_trans = exprdf_trans[, 2:ncol(exprdf_trans)]

ndnf_markers = aibs_primary_trans %>% mutate(ndnf = rowMeans(cbind(`Ndnf.Car4`, `Ndnf.Cxcl14`)), 
                                       l23 = rowMeans(cbind(`L2.Ngb`,`L2.3.Ptgs2`)),
                                       ratio = ndnf / l23) %>% 
  filter(ratio > 50, ndnf > 200) %>% arrange(-ndnf) %>% select(gene) %>% unlist() %>% as.character

sncg_markers = aibs_primary_trans %>% mutate(sncg = rowMeans(cbind(`Sncg`)), 
                                       l23 = rowMeans(cbind(`L2.Ngb`,`L2.3.Ptgs2`)),
                                       ratio = sncg / l23) %>% 
  filter(ratio > 50, sncg > 200) %>% arrange(-sncg) %>% select(gene) %>% unlist() %>% as.character

# filter out sparsely 
nnz_gene_counts = aibsExprDataDf %>% select(one_of(c(sncg_markers, 'primary_type'))) %>% filter(primary_type == 'Sncg') %>% select(-primary_type) %>% summarize_all(funs(sum(. > 1)))
sncg_markers = sncg_markers[which(nnz_gene_counts > 2)]

save(aibs_primary, ndnf_markers, sncg_markers, file = 'data/tasic_extra_markers.rda')



