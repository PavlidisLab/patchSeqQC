library(ggplot2)
library('biomaRt')
library(plyr)
library(dplyr)
library(tidyr)
library(R.matlab)
library(stringr)
library(cowplot)
library(ogbox)
library(readxl)
library(devtools)


#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=mart)

neuroelectro_ephys_vars = c("rin", "rmp","apthr","apamp",
                  "aphw", "tau", "ahpamp", "rheo","maxfreq", "cap", "adratio")



### LOAD Cadwell dataset

print('Loading Cadwell dataset')

cadwell_meta <- read.csv("data-raw/cadwell/E-MTAB-4092.sdrf.txt", sep = '\t')
cadwell_meta = cadwell_meta %>% separate(Source.Name, into = c('junk', 'num_idx'), sep = '_')

cadwell_ephys<- readMat("data-raw/cadwell/IntPhysiologyFinal.mat") %>% as.data.frame %>% t %>% as.data.frame
e = apply(cadwell_ephys, 2, function(x) as.numeric(unlist(x))) %>% as.data.frame
e$class = as.factor(cadwell_ephys$class %>% unlist)
e$typescore = as.factor(cadwell_ephys$typescore %>% unlist)

cadwell_ephys<- e

cadwell_ephys = cadwell_ephys %>% dplyr::rename(rmp = vrest, rin = inputres, apthr = thresh, aphw = width, apamp = amp, ahpamp = ahp)
consensus_ephys_props = neuroelectro_ephys_vars[neuroelectro_ephys_vars %in% colnames(cadwell_ephys)]

cadwell_ephys$num_idx = as.character(cadwell_ephys$idx)

cadwell_joined = dplyr::left_join(cadwell_meta, cadwell_ephys)
# cadwell_joined = merge(cadwell_meta, cadwell_ephys)

cadwell_expr = read.csv("data-raw/cadwell/tpmMatrix.genes", sep = '\t')
cadwell_expr = merge(annot, cadwell_expr, by.x = 'ensembl_gene_id', by.y = 'genes')
rownames(cadwell_expr) = make.names(cadwell_expr$mgi_symbol, unique = T)
cadwell_gene_names = rownames(cadwell_expr)

sample_inds = str_detect(colnames(cadwell_expr), 'ERR')
cadwell_expr = cadwell_expr[, sample_inds]

cadwell_expr = cadwell_expr %>% t() %>% as.data.frame

cadwell_expr = cadwell_expr + 1

cadwell_expr$Comment.ENA_RUN. = str_extract(rownames(cadwell_expr), 'ERR[0-9]+')
cadwell_expr$geo_id = cadwell_expr$Comment.ENA_RUN.
# cadwell_expr$Comment.ENA_RUN. = rownames(cadwell_expr)


# load count matrices to calculate ERCC percentages

cadwell_ercc = read.csv("data-raw/cadwell/countMatrix.ERCC", sep = '\t')
cadwell_count_matrix = read.csv("data-raw/cadwell/countMatrix.genes", sep = '\t')
cadwell_readcounts = read.csv("data-raw/cadwell/E-MTAB-4092.readcount", sep = ',', header = F)
cadwell_readcounts$geo_id = str_extract(cadwell_readcounts[, 1],'ERR[0-9]+')
colnames(cadwell_readcounts)[2] = 'read_count'

cadwell_readcounts = cadwell_readcounts %>% dplyr::select(geo_id, read_count)

ercc_sum = colSums(cadwell_ercc[, -1])
num_genes_count_matrix = colSums(cadwell_count_matrix[, -1] > 0)
read_count_mapped = colSums(cadwell_count_matrix[, -1])
ercc_df = cbind(ercc_sum, num_genes_count_matrix, read_count_mapped) %>% data.frame() %>% tibble::rownames_to_column(var = "geo_id")

ercc_df = merge(ercc_df, cadwell_readcounts)
ercc_df$ercc_pct = 100 * ercc_df$ercc_sum / ercc_df$read_count
ercc_df$ercc_pct_mapped = 100 * ercc_df$ercc_sum / ercc_df$read_count_mapped

ercc_df = ercc_df %>% dplyr::select(geo_id, ercc_sum, read_count, read_count_mapped, num_genes_count_matrix, ercc_pct, ercc_pct_mapped)

cadwell_expr = merge(cadwell_expr, ercc_df)


cadwell_joined = merge(cadwell_joined, cadwell_expr, by = 'Comment.ENA_RUN.')
rownames(cadwell_joined) = cadwell_joined$num_idx

cadwell_top_disc_genes = read.csv('~/ephys_analysis/data/cadwell_patchseq/top_disc_genes.csv')

# cadwell_joined[, 'colors'] = 'firebrick'
cadwell_joined[cadwell_joined$typescore %in% c(1, 2, 3), 'major_type'] = 'eNGC'
cadwell_joined[cadwell_joined$typescore %in% c(4, 5), 'major_type'] = 'SBC'
cadwell_joined[cadwell_joined$major_type %in% c('eNGC'), 'colors'] = '#F08080'
cadwell_joined[cadwell_joined$major_type %in% c('SBC'), 'colors'] = '#660066'

cadwell_joined$major_type = factor(cadwell_joined$major_type, levels = c('eNGC', 'Exc', 'Inh', 'SBC'))

cadwell_joined[cadwell_joined$Characteristics.cell.type. == 'Layer I Interneuron' & is.na(cadwell_joined$major_type), 'major_type'] = 'Inh'
cadwell_joined[cadwell_joined$Characteristics.cell.type. == 'Layer 2/3 Pyramidal' & is.na(cadwell_joined$major_type), 'major_type'] = 'Exc'

cadwell_joined$cell_id = make.names(cadwell_joined$major_type, unique = T)
cadwell_joined$sample_id = cadwell_joined$Comment.ENA_RUN.

cadwell_joined = cadwell_joined[!is.na(cadwell_joined$major_type), ] # remove 1 astrocyte sample

cadwell_joined = cadwell_joined[order(cadwell_joined$cell_id), ]

cadwell_ob = list(joined_df = cadwell_joined, gene_names = cadwell_gene_names, ephys_names = consensus_ephys_props)

### export data

cadwell_expr = cadwell_joined[, cadwell_gene_names] %>% t()
rownames(cadwell_expr) = make.names(cadwell_gene_names)
colnames(cadwell_expr) = cadwell_joined$geo_id
cadwell_expr = cadwell_expr -1

devtools::use_data(cadwell_expr, overwrite = T)

cadwell_meta = cadwell_joined[, c('geo_id', 'major_type', 'colors')]
cadwell_meta = dplyr::left_join(cadwell_meta, ercc_df %>% dplyr::select(geo_id, read_count, ercc_sum, ercc_pct)) %>% 
  dplyr::rename(sample_id = geo_id, ercc_count = ercc_sum)
rownames(cadwell_meta) = c()
devtools::use_data(cadwell_meta, overwrite = T)

### LOAD Foldy dataset
print('Loading Foldy dataset')


foldy_expr <- read.csv("data-raw/foldy/tpmMatrix.genes", sep = '\t')
soft_file = "data-raw/foldy/GSE75386_family.soft"

foldy_soft = softParser(soft_file)
foldy_soft = foldy_soft %>% dplyr::rename(sample_title = '!Sample_title', geo_id ='!Sample_geo_accession', sample_desc_name = '!Sample_source_name_ch1')


foldy_expr = merge(annot, foldy_expr, by.x = 'ensembl_gene_id', by.y = 'genes')
rownames(foldy_expr) = make.names(foldy_expr$mgi_symbol, unique = T)
foldy_gene_names = rownames(foldy_expr)
gsm_names = colnames(foldy_expr)[str_detect(colnames(foldy_expr), 'GSM')]
foldy_expr = foldy_expr[, gsm_names]
foldy_expr = foldy_expr %>% t() %>% as.data.frame

foldy_expr = foldy_expr + 1

foldy_expr$geo_id = rownames(foldy_expr)

foldy_expr = merge(foldy_soft %>% select(sample_title, geo_id, sample_desc_name) %>% dplyr::rename(cell_id = sample_title), foldy_expr, by = 'geo_id')
rownames(foldy_expr) = foldy_expr$geo_id

### load count matrices for ercc calculation


foldy_ercc <- read.csv("data-raw/foldy/countMatrix.ERCC", sep = '\t')
foldy_counts <- read.csv("data-raw/foldy/GSE75386_counts.genes", sep = '\t')

foldy_counts_combined = bind_rows(foldy_counts, foldy_ercc)
rownames(foldy_counts_combined) = make.names(foldy_counts_combined$genes, unique = T)

ercc_genes = rownames(foldy_counts_combined)[grepl('ERCC', rownames(foldy_counts_combined))] %>% as.character

foldy_readcounts = read.csv("data-raw/foldy/GSE75386.readcount", sep = ',', header = F)
foldy_readcounts$geo_id = str_extract(foldy_readcounts[, 1],'GSM[0-9]+')
colnames(foldy_readcounts)[2] = 'read_count'

foldy_readcounts = foldy_readcounts %>% select(geo_id, read_count) %>% distinct(geo_id, read_count)

ercc_sum = colSums(foldy_counts_combined[ercc_genes, -1])
ercc_sum[ercc_sum < 1000] = NA # if less than 1000 ercc counts detected, call NA
num_genes_count_matrix = colSums(foldy_counts[, -1] > 0)
ercc_df = cbind(ercc_sum, num_genes_count_matrix) %>% data.frame() %>% tibble::rownames_to_column(var = "geo_id")

ercc_df = merge(ercc_df, foldy_readcounts) 
ercc_df$ercc_pct = 100 * ercc_df$ercc_sum / ercc_df$read_count

ercc_df = ercc_df %>% select(geo_id, ercc_sum, read_count, ercc_pct)

foldy_expr = merge(foldy_expr, ercc_df)

# read in custom spreadsheet for ephys data, with cell id's provided by foldy's email
foldy_ephys_data <- read.delim("data-raw/foldy/foldy_ephys_data.csv")

foldy_ephys_data = foldy_ephys_data %>% separate(cell_id, c("num_id", "cell_id"), sep = " ")

foldy_ephys_data$rin = foldy_ephys_data$rin * 1000
foldy_ephys_data$aphw = foldy_ephys_data$aphw * 100
foldy_ephys_data$apwidth = foldy_ephys_data$apwidth * 100

consensus_ephys_props = neuroelectro_ephys_vars[neuroelectro_ephys_vars %in% colnames(foldy_ephys_data)]


# join ephys and expression data frames using cell id's as identifiers

foldy_j = right_join(foldy_ephys_data, foldy_expr, by= 'cell_id')
foldy_j$orig_cell_id = foldy_j$cell_id

foldy_j$major_type = factor(foldy_j$major_type, levels = c(foldy_j$major_type %>% levels, 'RS-PYR', 'BS-PYR'))

foldy_j[foldy_j$sample_desc_name == 'CA1 fast-spiking interneuron', 'major_type'] = 'FS-INT cells'
foldy_j[foldy_j$sample_desc_name == 'CA1 regular-spiking interneuron', 'major_type'] = 'RS-INT cells'
foldy_j[foldy_j$sample_desc_name == 'CA1 pyramidal cell', 'major_type'] = 'CA1-PYR cells'
foldy_j[foldy_j$sample_desc_name == 'Subiculum burst-spiking pyramidal cell', 'major_type'] = 'BS-PYR'
foldy_j[foldy_j$sample_desc_name == 'Subiculum regular-spiking pyramidal cell', 'major_type'] = 'RS-PYR'

foldy_j[foldy_j$major_type %in% c('FS-INT cells'), 'colors'] = 'red'
foldy_j[foldy_j$major_type %in% c('RS-INT cells'), 'colors'] = 'purple'
foldy_j[foldy_j$major_type %in% c('CA1-PYR cells'), 'colors'] = 'turquoise'
foldy_j[foldy_j$major_type %in% c('RS-PYR', 'BS-PYR'), 'colors'] = 'green'

foldy_j = foldy_j[!is.na(foldy_j$major_type), ] # remove bulk tissue samples

foldy_j$major_type = mapvalues(foldy_j$major_type, from = c("FS-INT cells", "RS-INT cells", "CA1-PYR cells"), to = c("FS-INT", "RS-INT", "CA1-PYR"))
foldy_j$major_type = factor(foldy_j$major_type, levels = c('BS-PYR', 'CA1-PYR', 'FS-INT', 'RS-INT', 'RS-PYR'))

foldy_joined = foldy_j

foldy_joined$cell_id = make.names(foldy_joined$major_type, unique = T)
foldy_joined$sample_id = foldy_joined$geo_id

foldy_joined = foldy_joined[order(foldy_joined$cell_id), ]

foldy_ob = list(joined_df = foldy_joined, gene_names = foldy_gene_names, ephys_names = consensus_ephys_props)


### Load Fuzik

print('Loading Fuzik dataset')

soft_file = 'data-raw/fuzik/GSE70844_family.soft'
fuzik_meta = softParser(soft_file)

cell_id_names = fuzik_meta$`!Sample_characteristics_ch1 = cell id`

ephys_file = 'data-raw/fuzik/fuzik_ephys.csv'
fuzik_ephys = read.csv(ephys_file)

fuzik_ephys_cleaned = fuzik_ephys

fuzik_quant_ephys_vals = colnames(fuzik_ephys)[4:96]
fuzik_ephys_cleaned[, 4:96] <- lapply(fuzik_ephys[, fuzik_quant_ephys_vals], function(x) as.numeric(as.character(x)))
fuzik_ephys_cleaned[is.na(fuzik_ephys_cleaned)] = NaN


fuzik_ephys_cleaned$updated_ids = sapply(fuzik_ephys_cleaned$cell_id, str_replace, 'C', 'cell') %>% tolower %>% str_trim

new_ids = sapply(fuzik_ephys_cleaned$cell_id, str_replace, 'C', 'cell') %>% tolower %>% str_trim

matching_inds = sapply(new_ids, FUN = function(x1) match(x1, str_trim(tolower(cell_id_names))))

rownames(fuzik_meta)[matching_inds]

fuzik_ephys_cleaned$expr_ids = fuzik_meta$`!Sample_title`[matching_inds] %>% str_replace('_single cell', '')
fuzik_ephys_cleaned$GSM_id = rownames(fuzik_meta)[matching_inds]

fuzik_ephys_cleaned = fuzik_ephys_cleaned %>% dplyr::rename(rmp = Vrest, rin = R.input..MOhm., tau = Tau..ms.,
                                                            maxfreq = max.Frequency..Hz., apthr = AP.threshold..mV., 
                                                            aphw = AP.halfwidth..ms., apamp = AP.ampl..mV., ahpamp = AHP.ampl..mV., 
                                                            rheo = Rheobasic.current.step, adratio = Adaptation.2ndlast.2XTHR) %>%
  mutate(cap = tau / rin * 1000 , adratio = 1/adratio)
fuzik_ephys_cleaned$adratio[is.infinite(fuzik_ephys_cleaned$adratio)] = NA # map inf adratios to NA

consensus_ephys_props = neuroelectro_ephys_vars

# fuzik_ephys_cleaned %>% select()


# predict maxfreq values which are 0 based on Avg.frequency.2X.THR..Hz.

# fit lm to model maxfreq using Avg Freq
maxfreq_model = lm(maxfreq ~ Avg.frequency.2X.THR..Hz. - 1, data = fuzik_ephys_cleaned %>% filter(maxfreq > 0, Avg.frequency.2X.THR..Hz. > 0))

# update 0 maxfreq inds using model predictions
update_inds = fuzik_ephys_cleaned$maxfreq == 0 & !is.na(fuzik_ephys_cleaned$maxfreq)
fuzik_ephys_cleaned[update_inds, 'maxfreq'] = predict(maxfreq_model, fuzik_ephys_cleaned[update_inds, ])

# since since some maxfreq inds are still 0 call these NaN cause they can't be 0
update_inds = fuzik_ephys_cleaned$maxfreq == 0 & !is.na(fuzik_ephys_cleaned$maxfreq)
fuzik_ephys_cleaned[update_inds, 'maxfreq'] = NaN

# since since some maxfreq inds are still 0 call these NaN cause they can't be 0
update_inds = fuzik_ephys_cleaned$Avg.frequency.2X.THR..Hz. == 0 & !is.na(fuzik_ephys_cleaned$Avg.frequency.2X.THR..Hz.)

fuzik_ephys_cleaned[update_inds, 'Avg.frequency.2X.THR..Hz.'] = NaN

# load expression matrix and transpose 

expr_file = 'data-raw/fuzik/GSE70844_Fuzik_et_al_molcounts.xls'
fuzik_expr = read_excel(expr_file)
colnames(fuzik_expr)[1] = 'gene'

fuzik_expr_t = fuzik_expr[, -1] %>% t() %>% as.data.frame()

# calculate total molecules detected per sample
cell_total_counts = fuzik_expr_t %>% rowSums()

# normalize by total counts to molecules per million
fuzik_expr_t = fuzik_expr_t *1E6 / cell_total_counts
fuzik_expr_t = fuzik_expr_t + 1

colnames(fuzik_expr_t) = lapply(fuzik_expr[, 1], make.names) %>% unlist
fuzik_gene_names = colnames(fuzik_expr_t)


fuzik_expr_t$sample.name = rownames(fuzik_expr_t)

fuzik_joined = merge(fuzik_expr_t, fuzik_ephys_cleaned, by.x = 'sample.name', by.y = 'expr_ids')
fuzik_joined$sub_type = fuzik_joined$Cell.type

fuzik_joined[str_detect(fuzik_joined$sub_type, 'L'), 'colors'] = 'turquoise'
fuzik_joined[str_detect(fuzik_joined$sub_type, 'T'), 'colors'] = 'purple'

fuzik_joined[str_detect(fuzik_joined$sub_type, 'L'), 'major_type'] = 'Exc'
fuzik_joined[str_detect(fuzik_joined$sub_type, 'T'), 'major_type'] = 'Inh'
fuzik_joined$major_type = factor(fuzik_joined$major_type, levels = c('Exc', 'Inh'))
# fuzik_joined$broad_type = df <- data.frame(matrix(unlist(fuzik_joined[, consensus_ephys_props]), nrow=nrow(fuzik_joined), byrow=F), stringsAsFactors = F)
# fuzik_joined[, consensus_ephys_props] = df
fuzik_joined = fuzik_joined %>% arrange(major_type) # sort by cell type

fuzik_joined$cell_id = make.names(fuzik_joined$major_type, unique = T)
fuzik_joined$sample_id = fuzik_joined$sample.name

fuzik_joined = fuzik_joined[order(fuzik_joined$cell_id), ]

fuzik_ob = list(joined_df = fuzik_joined, gene_names = fuzik_gene_names, ephys_names = consensus_ephys_props)

patch_seq_datasets = list(cadwell = cadwell_ob, foldy = foldy_ob, fuzik = fuzik_ob)

print('saving all patch seq datasets to an R object, patch_seq_datasets')
saveRDS(patch_seq_datasets, file = 'data/patch_seq_datasets.rda')
