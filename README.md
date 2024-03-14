# CHAI: consensus Clustering tHrough similArIty matrIx integratIon for single cell type identification

### Introduction 
CHAI (consensus Clustering tHrough similArIty matrIces for single cell type identification) is a consensus clustering framework that offers two methods for consensus clustering: Average Similarity (AvgSim) and Similarity Network Fusion (SNF).

![vr_flowchart](https://github.com/lodimk2/chai/assets/69815640/21202365-38f7-4fa9-aeff-c8f5b14c9fe9)

### Installation Instructions 

You may install CHAI using devtools. 

```devtools::install_github("lodimk2/chai")```

### Quick Start

CHAI contains two methods for consensus clustering: CHAI-AvgSim and CHAI-SNF. We provide wrapper functions to run either of these. 

We provide example data from the Baron Mouse 1 Dataset (Baron et al., 2016). Data should be in the form of a SingleCellExperiment object, with counts and logcounts defined. 

#### Load Dependencies:

```
library(chai)
library(scSHC)
library(RaceID)
library(SC3)
library(SingleCellExperiment)
library(CHOIR)
```

#### Load Data:

```
data("baron_mouse_1")
# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(baron_mouse_1)))
# Add logcounts 
sce <- scuttle::logNormCounts(sce)
```
#### CHAI-AvgSim:

```
# If "eval" is set to TRUE, CHAI will evaluate the best_k for Spectral Clustering using silhouette score. If "eval" is set to false, CHAI will use best_k as the k value for silhouette score.
# The default is eval being set to TRUE
best_k <- 15
CHAI_AvgSim <- function(sce,best_k,eval = TRUE)
```

#### CHAI-SNF:

```
# If "eval" is set to TRUE, CHAI will evaluate the best_k for Spectral Clustering using silhouette score. If "eval" is set to false, CHAI will use best_k as the k value for silhouette score.
# The default is eval being set to TRUE
best_k <- 15
CHAI_SNF <- function(sce,best_k,eval = TRUE)
```
### In Depth Tutorial
For a more detailed tutorial, including adding your own clustering algorithms to CHAI and including other "omics" data, please see the vignette. 

### References:
1. Baron M, Veres A, Wolock SL, Faust AL, Gaujoux R, Vetere A, Ryu JH, Wagner BK, Shen-Orr SS, Klein AM, Melton DA, Yanai I. A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell Syst. 2016 Oct 26;3(4):346-360.e4. doi: 10.1016/j.cels.2016.08.011. Epub 2016 Sep 22. PMID: 27667365; PMCID: PMC5228327. 
