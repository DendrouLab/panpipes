[![PyPI version](https://badge.fury.io/py/panpipes.svg)](https://badge.fury.io/py/panpipes)

# Panpipes - multimodal single cell pipelines

## Overview

Panpipes is a set of computational workflows designed to automate multimodal single-cell and spatial transcriptomic analyses by incorporating widely-used Python-based tools to perform quality control, preprocessing, integration, clustering, and reference mapping at scale. 
Panpipes allows reliable and customisable analysis and evaluation of individual and integrated modalities, thereby empowering decision-making before downstream investigations.

**See our [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/) and our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.11.532085v2)**:  

These workflows make use of [cgat-core](https://github.com/cgat-developers/cgat-core):

Available workflows:

1. "ingest" : Ingest data and compute quality control metrics
2. "preprocess" : Filter and normalize per modality
3. "integration" : Integrate and batch correct using single and multimodal methods
4. "clustering" : Cluster cells per modality
5. "refmap" : Map queries against reference datasets
6. "vis" : Visualize metrics from other pipelines in the context of experiment metadata
7. "qc_spatial" : Ingest spatial transcriptomics data (Vizgen, Visium) and compute quality control metrics
8. "preprocess_spatial" : Filtering and normalize spatial transcriptomics data
9. "deconvolution_spatial" : Deconvolve cell types of spatial transcriptomics slides

## Installation and configuration

See [installation instructions here](https://panpipes-pipelines.readthedocs.io/en/latest/install.html)

Oxford BMRC Rescomp users find additional advice in [docs/installation_rescomp](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md)

## Releases

`panpipes v0.4.0` is out [now](./CHANGELOG.md)!

The `ingest` workflow now expects different headers for the RNA and Protein modalities.
Check the example [submission file](https://github.com/DendrouLab/panpipes/blob/main/docs/usage/sample_file_qc_mm.md) and the [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/usage/setup_for_qc_mm.html) for more detailed instructions.

## Citation

[Panpipes: a pipeline for multiomic single-cell and spatial transcriptomic data analysis
Fabiola Curion, Charlotte Rich-Griffin, Devika Agarwal, Sarah Ouologuem, Tom Thomas, Fabian J. Theis, Calliope A. Dendrou
bioRxiv 2023.03.11.532085; doi: https://doi.org/10.1101/2023.03.11.532085](https://www.biorxiv.org/content/10.1101/2023.03.11.532085v2)

## Contributors

Created and Maintained by Charlotte Rich-Griffin and Fabiola Curion.
Additional contributors: Sarah Ouologuem, Devika Agarwal, and Tom Thomas.
