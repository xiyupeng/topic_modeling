Latent Dirichlet Allocation with application on Longitudinal Flow Cytometry Data
=======

This repository contains all codes used in the paper: 
“Uncovering the hidden structure of dynamic T cell
composition in peripheral blood during cancer immunotherapy: a topic modeling
approach”

# Prerequisites

Below are three core R packages we used for LDA analysis

```
library(topicmodels)
library(slam)
library(tidytext)
```

We use the R package *topicmodels* for model inference, R package *slam* for preparing the required input data, R package *tidytext* for extracting the output data. 
The three R packages are important for reproducing the LDA analysis we described in the paper on other datasets. 

There are also other R packages we used in the script, like *Seurat* package for UMAP visualization and clustering, *ComplexHeatmap* for heatmaps. 
you can also use other clustering methods, like FlowSOM, or visualization tools instead.


# Preparing input

In text mining, LDA requires the input as a document-by-term count matrix, where each row represent each document,each column represent each term, each entry in the matrix is the number of occurrences of each term (a word is a single occurrence of a term). Motivated by the similarities between text data mining and single-cell analysis, for single-cell analysis, LDA consider cells as words, cell types as terms, patient samples as documents, biological processes as topics.

Before applying the LDA model to single-cell datasets, we need to prepare the cell type-by-sample count matrix as the input of LDA. 
One common approach to obtain the cell type count matrix is to pool all cells together and do the clustering (Before doing the pooled clustering, you may want to check if there is a batch effect between samples). 
`X50_single_cell_analysis.R` contains all codes that we used to cluster the 17M+ cells.
The Louvain method in *Seurat* package was used to cluster cells and prepare the cell type-by-sample matrix. 
Cell types were manually annotated based on their marker expression.

# Usage

For more details of the usage of the LDA method,
please check the [tutorial](https://xiyupeng.github.io/LDA_examples/), where we show two examples of the 
application of LDA on single-cell dataset:

- scRNA-seq data of liver cancer patients. (Data from a Nature [paper](https://www.nature.com/articles/s41586-022-05400-x#Bib1))
- Longitudinal flow cytometry data of melanoma patients. (the paper)

# Data

Additional data used to generate figures were provided in:

Peng, Xiyu (2023), “flow cytometry dataset of melanoma patients”, Mendeley Data, V1, doi: 10.17632/d7nkgfhc8z.1

# Citation

X. Peng, J. Lee, M. Adamow, C. Maher, M. A. Postow, M. Callahan, K. S.
Panageas, R. Shen (2022+). “Uncovering the hidden structure of dynamic T cell
composition in peripheral blood during cancer immunotherapy: a topic modeling
approach”. (Under Review)

# Contact us

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com, pengx1@mskcc.org)

