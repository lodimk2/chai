# CHAI: consensus Clustering tHrough similArIty matrix integratIon for single cell type identification 
<div style="text-align: center;">
    <img src="https://github.com/lodimk2/chai/assets/69815640/c3894c8e-91c3-45a5-acaa-0a6ca8eb79e6" alt="vr_chai_logi" width="250px" height="300px">
</div>


## Latest News
### Version 0.99.0 (2023-03-15)
- Initial release of the package.
- Submitted to Bioconductor.

For the full log of news and updates, please check the [NEWS.md](NEWS.md) file.

## Introduction 
CHAI (consensus Clustering tHrough similArIty matrix integratIon for single cell type identification) is a consensus clustering framework that offers two methods for consensus clustering: Average Similarity (AvgSim) and Similarity Network Fusion (SNF) (Wang et al. 2014).

![vr_flowchart](https://github.com/lodimk2/chai/assets/69815640/21202365-38f7-4fa9-aeff-c8f5b14c9fe9)

## Installation Instructions 

You may install CHAI using devtools. 

```devtools::install_github("lodimk2/chai")```

## Quick Start

CHAI contains two methods for consensus clustering: CHAI-AvgSim and CHAI-SNF. We provide wrapper functions to run either of these. 

We provide example data from the Baron Mouse 1 Dataset (Baron et al., 2016). Data should be in the form of a SingleCellExperiment object, with counts and logcounts defined. 

### Load Dependencies:

```
library(chai)
library(scSHC)
library(RaceID)
library(SC3)
library(SingleCellExperiment)
library(CHOIR)
```

### Load Data:

```
data("baron_mouse_1")
# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(baron_mouse_1)))
# Add logcounts 
sce <- scuttle::logNormCounts(sce)
```
### CHAI-AvgSim:

```
# If "eval" is set to TRUE, CHAI will evaluate the best_k for Spectral Clustering using silhouette score. If "eval" is set to false, CHAI will use best_k as the k value for silhouette score.
# The default is eval being set to TRUE
best_k <- 15
sce <- CHAI_AvgSim <- function(sce,best_k,eval = TRUE)
```

### CHAI-SNF:

```
# If "eval" is set to TRUE, CHAI will evaluate the best_k for Spectral Clustering using silhouette score. If "eval" is set to false, CHAI will use best_k as the k value for silhouette score.
# The default is eval being set to TRUE
best_k <- 15
sce <- CHAI_SNF <- function(sce,best_k,eval = TRUE)
```
## In Depth Tutorial
For a more detailed tutorial, including adding your own clustering algorithms to CHAI and including other "omics" data, please see the vignette [inst/chai.html](inst/chai.html). 


## References:
1. Baron M, Veres A, Wolock SL, Faust AL, Gaujoux R, Vetere A, Ryu JH, Wagner BK, Shen-Orr SS, Klein AM, Melton DA, Yanai I. A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell Syst. 2016 Oct 26;3(4):346-360.e4. doi: 10.1016/j.cels.2016.08.011. Epub 2016 Sep 22. PMID: 27667365; PMCID: PMC5228327.
2. Wang, B., Mezlini, A., Demir, F. et al. Similarity network fusion for aggregating data types on a genomic scale. Nat Methods 11, 333–337 (2014). https://doi.org/10.1038/nmeth.2810

## Citation:
If you use CHAI, please cite our manuscript: 

Lodi, M., Lodi, M., Osei, K., Ranganathan, V., Hwang, P., & Ghosh, P. (2024). CHAI: Consensus Clustering Through Similarity Matrix Integration for Cell-Type Identification. bioRxiv. https://doi.org/10.1101/2024.03.19.585758

```
@article {Lodi2024.03.19.585758,
	author = {Musaddiq Lodi and Muzammil Lodi and Kezie Osei and Vaishnavi Ranganathan and Priscilla Hwang and Preetam Ghosh},
	title = {CHAI: Consensus Clustering Through Similarity Matrix Integration for Cell-Type Identification},
	elocation-id = {2024.03.19.585758},
	year = {2024},
	doi = {10.1101/2024.03.19.585758},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Several methods have been developed to computationally predict cell-types for single cell RNA sequencing (scRNAseq) data. As methods are developed, a common problem for investigators has been identifying the best method they should apply to their specific use-case. To address this challenge, we present CHAI (consensus Clustering tHrough similArIty matrix integratIon for single cell type identification), a wisdom of crowds approach for scRNAseq clustering. CHAI presents two competing methods which aggregate the clustering results from seven state of the art clustering methods: CHAI-AvgSim and CHAI-SNF. Both methods demonstrate improved performance on a diverse selection of benchmarking datasets, besides also outperforming a previous consensus clustering method. We demonstrate CHAI{\textquoteright}s practical use case by identifying a leader tumor cell cluster enriched with CDH3. CHAI provides a platform for multiomic integration, and we demonstrate CHAI-SNF to have improved performance when including spatial transcriptomics data. CHAI is intuitive and easily customizable; it provides a way for users to add their own clustering methods to the pipeline, or down-select just the ones they want to use for the clustering aggregation. CHAI is available as an open source R package on GitHub: https://github.com/lodimk2/chaiCompeting Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2024/03/22/2024.03.19.585758},
	eprint = {https://www.biorxiv.org/content/early/2024/03/22/2024.03.19.585758.full.pdf},
	journal = {bioRxiv}
}
```

## Questions and Support:
For any questions or comments, please reach out to Musaddiq Lodi @ lodimk2@vcu.edu

