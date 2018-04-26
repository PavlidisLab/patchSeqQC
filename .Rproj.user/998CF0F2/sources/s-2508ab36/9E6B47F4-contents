### master file for running all scripts and generating all figures

# big make file for patch seq qc paper

### load patch seq datasets from csv files
source('analysis/LoadingScripts/LoadPatchSeqData.R')

### load dissociated cell reference data from .rda files
source('analysis/LoadingScripts/LoadDissociatedCellData.R')
# 
# source('analysis/LoadingScripts/LoadTasic2016.R')
# source('analysis/LoadingScripts/LoadZeisel.R')

### map cell types across datasets
source('analysis/mapCellTypes.R')

### load marker genes and update markers using expression in Tasic 2016 and Zeisel 2015 data
source('analysis/defineMarkers.R')

### write mappings between broad cell types
source('analysis/writeCellTypeMappings.R')

### calculate marker expression for all dissociated cells and patch seq samples
source('analysis/calculateMarkerExpression.R') 

### plot heatmaps and generate figures needed for figure 1, illustrating contamination issue in patch seq data
source('analysis/illustrateContamination.R')

### count numbers of genes expressed in patch seq data and how this compares with ercc ratio and contamination
source('analysis/numGenesPerDataset.R')

### calculate corrlations between ephys and expression for tasic/aibs cell types data
# source('analysis/LoadingScripts/LoadAibsEphys.R')
# source('analysis/mergeAibsExprEphys.R') 
source('analysis/calcAibsSigGenes.R')

### calculate correlations between ephys and expression data for patch-seq data
source('analysis/ephysComparisionAnalysis.R')

### plot final ephys expression comparison fig
source('analysis/ephysComparisionAnalysisAgg.R')
