# load data from Zeisel dissociated cell brain dataset

# largely copied from ogan's version of the same function: https://github.com/PavlidisLab/neuroExpressoAnalysis/blob/master/analysis/00.DownloadPreprocess/02.BrainSingleCellData.R

library(XLConnect)
library(magrittr)
library(ogbox)
library(dplyr)

print('Loading Zeisel 2015 data')

# # mouse rna seq data download Zeisel et al. download------
# dir.create('data-raw/ZeiselMouse', showWarnings=FALSE)
# download.file(url= 'http://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt', 
#               destfile='data-raw/ZeiselMouse/mouseRNASeq_Zeisel 2015.txt')
# download.file(url = 'http://science.sciencemag.org/highwire/filestream/628248/field_highwire_adjunct_files/1/aaa1934_TableS1.xlsx',
#               destfile = 'data-raw/ZeiselMouse/markerGenes.xlsx')
# download.file(url = 'http://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt',
#               destfile = 'data-raw/ZeiselMouse/expression_spikes_17-Aug-2014.txt')
# 

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}


# mouse rna seq data download Zeisel et al. process------
rnaSeq = read.table('data-raw/ZeiselMouse/mouseRNASeq_Zeisel 2015.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
rnaMeta = rnaSeq[1:10,3:ncol(rnaSeq)]
rnaMeta = as.data.frame(t(rnaMeta))
colnames(rnaMeta) = rnaSeq[1:10,2]

rnaExp = rnaSeq[12:nrow(rnaSeq),3:ncol(rnaSeq)]
rnaExp = apply(rnaExp,2,as.numeric)
rnaExp = matrix(unlist(rnaExp), nrow = nrow(rnaExp))
rownames(rnaExp) = rnaSeq[12:nrow(rnaSeq),1]
# rnaCelIDs = as.numeric(as.character(rnaSeq[12:nrow(rnaSeq), 2]))
# rnaExp = rnaExp[,rnaMeta$tissue %in% 'sscortex']
# rnaMeta = rnaMeta[rnaMeta$tissue %in% 'sscortex',]
# remove low expressed ones 
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp >= 1,]

rnaMeta$tissue %<>% as.character
rnaMeta$`group #` %<>% as.numeric.factor
rnaMeta$`total mRNA mol` %<>% as.numeric.factor
rnaMeta$well %<>% as.numeric.factor
rnaMeta$sex %<>% as.numeric.factor
rnaMeta$age %<>% as.numeric.factor
rnaMeta$diameter %<>% as.numeric.factor
rnaMeta$cell_id %<>% as.char
rnaMeta$level1class %<>% as.char
rnaMeta$level2class %<>% as.char

ZeiselMouseExp = rnaExp
ZeiselMouseMeta = rnaMeta

zeisel_cell_types = ZeiselMouseMeta$level2class %>% unique
ZeiselMouseExp_subgroup_mean = lapply(zeisel_cell_types, function(cell_type_name){
  cell_inds = which(ZeiselMouseMeta$level2class == cell_type_name)
  mean_exprs = rowMeans(ZeiselMouseExp[, cell_inds])
  return(mean_exprs)
})
names(ZeiselMouseExp_subgroup_mean) = zeisel_cell_types
zeisel_grps_mean = bind_cols(ZeiselMouseExp_subgroup_mean) %>% as.data.frame()
rownames(zeisel_grps_mean) = rownames(ZeiselMouseExp)

zeiselExprDataDf = cbind(ZeiselMouseMeta, t(ZeiselMouseExp))
zeisel_gene_names = rownames(ZeiselMouseExp)

saveRDS(zeisel_grps_mean, file = 'data/zeisel_group_mean.rda')

saveRDS(zeiselExprDataDf, file = 'data/zeisel_expression.rda')
saveRDS(zeisel_gene_names, file = 'data/zeisel_gene_names.rda')



