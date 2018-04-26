library(ggplot2)
library(dplyr)
library(tidyr)
library(R.matlab)
library(stringr)
source('~/ephys_analysis/ephysExprPlots.R')
library(cowplot)
library(ogbox)


#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=mart)



cadwell_meta <- read.csv("data-raw/cadwell/E-MTAB-4092.sdrf.txt", sep = '\t')
cadwell_meta = cadwell_meta %>% separate(Source.Name, into = c('junk', 'num_idx'), sep = '_')

cadwell_ephys<- readMat("data-raw/cadwell/IntPhysiologyFinal.mat") %>% as.data.frame %>% t %>% as.data.frame
e = apply(cadwell_ephys, 2, function(x) as.numeric(unlist(x))) %>% as.data.frame
e$class = as.factor(cadwell_ephys$class %>% unlist)
e$typescore = as.factor(cadwell_ephys$typescore %>% unlist)

cadwell_ephys<- e

cadwell_ephys = cadwell_ephys %>% dplyr::rename(rmp = vrest, rin = inputres, apthr = thresh, aphw = width, apamp = amp, ahpamp = ahp)
consensus_ephys_props = response_vars[response_vars %in% colnames(cadwell_ephys)]
consensus_ephys_props_cadwell = response_vars[response_vars %in% colnames(cadwell_ephys)]


cadwell_joined = merge(cadwell_meta, cadwell_ephys, by.x = 'num_idx', by.y = 'idx')


cadwell_expr = read.csv("data-raw/cadwell/tpmMatrix.genes", sep = '\t')
cadwell_expr = merge(annot, cadwell_expr, by.x = 'ensembl_gene_id', by.y = 'genes')
rownames(cadwell_expr) = make.names(cadwell_expr$mgi_symbol, unique = T)
cadwell_gene_names = rownames(cadwell_expr)

sample_inds = str_detect(colnames(cadwell_expr), 'ERR')
cadwell_expr = cadwell_expr[, sample_inds]

cadwell_expr = cadwell_expr %>% t() %>% as.data.frame

cadwell_expr = cadwell_expr + 1

cadwell_expr$Comment.ENA_RUN. = str_extract(rownames(cadwell_expr), 'ERR[0-9]+')
# cadwell_expr$Comment.ENA_RUN. = rownames(cadwell_expr)

cadwell_joined = merge(cadwell_joined, cadwell_expr, by = 'Comment.ENA_RUN.')
rownames(cadwell_joined) = cadwell_joined$num_idx

cadwell_top_disc_genes = read.csv('~/ephys_analysis/data/cadwell_patchseq/top_disc_genes.csv')

# cadwell_joined[, 'colors'] = 'firebrick'
cadwell_joined[cadwell_joined$typescore %in% c(1, 2, 3), 'major_type'] = 'eNGC'
cadwell_joined[cadwell_joined$typescore %in% c(4, 5), 'major_type'] = 'SBC'
cadwell_joined[cadwell_joined$major_type == 'eNGC', 'colors'] = '#F08080'
cadwell_joined[cadwell_joined$major_type == 'SBC', 'colors'] = '#660066'

cadwell_joined$cell_id = make.names(cadwell_joined$major_type, unique = T)
cadwell_joined$sample_id = cadwell_joined$Comment.ENA_RUN.

cadwell_joined = cadwell_joined[order(cadwell_joined$cell_id), ]