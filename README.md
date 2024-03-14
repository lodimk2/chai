# CHAI: consensus Clustering tHrough similArIty matrIx integratIon for single cell type identification

### Introduction 
CHAI (consensus Clustering tHrough similArIty matrIces for single cell type identification) is a consensus clustering framework that offers two methods for consensus clustering: Average Similarity (AvgSim) and Similarity Network Fusion (SNF).

![vr_flowchart](https://github.com/lodimk2/chai/assets/69815640/21202365-38f7-4fa9-aeff-c8f5b14c9fe9)

### Installation Instructions 

You may install CHAI using devtools. 

```devtools::install_github("lodimk2/chai")```

### Quick Start

CHAI contains two methods for consensus clustering: CHAI-AvgSim and CHAI-SNF. We provide wrapper functions to run either of these. 

We provide example data from the Baron Mouse 1 Dataset [citation needed here]. Data should be in the form of a SingleCellExperiment object, with counts and logcounts defined. 

#### Load Data:

```
data("baron_mouse_1")
# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(baron_mouse_1)))
# Add logcounts 
sce <- scuttle::logNormCounts(sce)
```


### Tutorial
For a tutorial, please see the vignettes directory. 
