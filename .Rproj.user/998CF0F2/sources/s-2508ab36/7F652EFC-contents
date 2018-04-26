
### annotate cell types in tasic and zeisel to standardized nomenclature based on neuroexpresso names


# assign and standardize major type names in tasic cells
tasic_major_types = aibsExprDataDf$primary_type %>% as.character()
tasic_major_types[str_detect(tasic_major_types, 'Micr')] = "Microglia"
tasic_major_types[str_detect(tasic_major_types, 'Astr')] = "Astrocyte"
tasic_major_types[str_detect(tasic_major_types, 'Oligo')] = "Oligodendrocyte"
tasic_major_types[str_detect(tasic_major_types, 'OPC')] = "OPC"
tasic_major_types[str_detect(tasic_major_types, 'Endo')] = "Endothelial"
tasic_major_types[str_detect(aibsExprDataDf$broad_type, 'Gluta')] = "Pyramidal"
tasic_major_types[str_detect(aibsExprDataDf$broad_type, 'GABA')] = "Inhibitory"
tasic_broad_types = factor(tasic_major_types)

tasic_major_types[str_detect(aibsExprDataDf$primary_type, 'Sncg')] = "Sncg"
tasic_major_types[str_detect(aibsExprDataDf$primary_type, 'Ndnf')] = "Ndnf"
tasic_major_types[str_detect(aibsExprDataDf$primary_type, 'Pvalb')] = "Pvalb"
tasic_sub_types = factor(tasic_major_types)

aibsExprDataDf$norm_sub_types = tasic_sub_types
aibsExprDataDf$norm_broad_types = tasic_broad_types



#### calculate expected marker sums in Zeisel data

zeisel_major_types = zeiselExprDataDf$level2class  %>% as.character()
zeisel_major_types[str_detect(zeisel_major_types, 'Mgl')] = "Microglia"
zeisel_major_types[str_detect(zeisel_major_types, 'Astr')] = "Astrocyte"
zeisel_major_types[str_detect(zeisel_major_types, 'Oligo')] = "Oligodendrocyte"
# zeisel_major_types[str_detect(zeisel_major_types, 'OPC')] = "OPC"
zeisel_major_types[str_detect(zeisel_major_types, 'Vsmc')] = "Endothelial"
zeisel_major_types[str_detect(zeisel_major_types, 'Pyr') & ! str_detect(zeisel_major_types, 'Int')] = "Pyramidal"
zeisel_major_types[str_detect(zeisel_major_types, 'Int') & ! str_detect(zeisel_major_types, 'Pyr')] = "Inhibitory"
zeisel_broad_types = factor(zeisel_major_types)


# map zeisel interneuron types onto tasic types based on mappings provided by MetaNeighbor
zeisel_major_types[zeiselExprDataDf$level2class %in% c('Int15', 'Int12')] = "Ndnf"
zeisel_major_types[str_detect(zeiselExprDataDf$level2class, 'Int5')] = "Sncg"
zeisel_major_types[str_detect(zeiselExprDataDf$level2class, 'Int3')] = "Pvalb"

zeisel_sub_types = factor(zeisel_major_types)

zeiselExprDataDf$norm_broad_type = zeisel_broad_types
zeiselExprDataDf$norm_sub_type = zeisel_sub_types
# 
# zeiselExprDataDf = cbind(ZeiselMouseMeta, t(ZeiselMouseExp))
# zeisel_gene_names = rownames(ZeiselMouseExp)


# map broad cell type (contam_type) to each defined_major type in patch seq datasets

patch_seq_datasets$cadwell$joined_df$contam_type = 'Ndnf'
patch_seq_datasets$cadwell$joined_df[str_detect(patch_seq_datasets$cadwell$joined_df$major_type, 'Exc'), 'contam_type'] = 'Pyramidal'
patch_seq_datasets$cadwell$joined_df$contam_type = factor(patch_seq_datasets$cadwell$joined_df$contam_type)

patch_seq_datasets$foldy$joined_df$contam_type = 'Sncg'
patch_seq_datasets$foldy$joined_df[str_detect(patch_seq_datasets$foldy$joined_df$major_type, 'RS-INT'), 'contam_type'] = 'Sncg'
patch_seq_datasets$foldy$joined_df[str_detect(patch_seq_datasets$foldy$joined_df$major_type, 'FS-INT'), 'contam_type'] = 'Pvalb'
patch_seq_datasets$foldy$joined_df[str_detect(patch_seq_datasets$foldy$joined_df$major_type, 'PYR'), 'contam_type'] = 'Pyramidal'

patch_seq_datasets$foldy$joined_df$contam_type = factor(patch_seq_datasets$foldy$joined_df$contam_type)

patch_seq_datasets$fuzik$joined_df$contam_type = 'Pyramidal'
patch_seq_datasets$fuzik$joined_df[str_detect(patch_seq_datasets$fuzik$joined_df$major_type, 'Exc'), 'contam_type'] = 'Pyramidal'
patch_seq_datasets$fuzik$joined_df[str_detect(patch_seq_datasets$fuzik$joined_df$major_type, 'Inh'), 'contam_type'] = 'Ndnf'

patch_seq_datasets$fuzik$joined_df$contam_type = factor(patch_seq_datasets$fuzik$joined_df$contam_type)


