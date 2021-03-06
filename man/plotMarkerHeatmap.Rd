% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMarkerHeatmaps.R
\name{plotMarkerHeatmap}
\alias{plotMarkerHeatmap}
\title{Generate a heatmap showing top n marker expression per cell type per patch seq dataset}
\usage{
plotMarkerHeatmap(markerlist, expr_matrix, show_legend = T,
  show_cell_labels = F)
}
\arguments{
\item{markerlist}{list of lists with marker genes per cell type}

\item{expr_matrix}{gene expression matrix (samples x genes)}

\item{show_legend}{boolean to show gene expression heatmap legend}

\item{show_cell_labels}{boolean to show cell labels}
}
\value{
a pheatmap object with resulting heatmap
}
\description{
Generate a heatmap showing top n marker expression per cell type per patch seq dataset
}
