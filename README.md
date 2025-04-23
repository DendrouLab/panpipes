[![PyPI version](https://badge.fury.io/py/panpipes.svg)](https://badge.fury.io/py/panpipes)

# Panpipes - multimodal single cell pipelines

## Overview

Panpipes is a set of computational workflows designed to automate multimodal single-cell and spatial transcriptomic analyses by incorporating widely-used Python-based tools to perform quality control, preprocessing, integration, clustering, and reference mapping at scale.
Panpipes allows reliable and customisable analysis and evaluation of individual and integrated modalities, thereby empowering decision-making before downstream investigations.

**See our [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/)**  

**Panpipes is on [Genome Biology!](https://link.springer.com/article/10.1186/s13059-024-03322-7)**

These workflows make use of [cgat-core](https://github.com/cgat-developers/cgat-core)

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

For detailed installation instructions (including those for Apple Silicon machines), refer to the [installation instructions here](https://panpipes-pipelines.readthedocs.io/en/latest/install.html).

We recommend installing panpipes in a conda environment.
We provide a minimal conda config file in `pipeline_env.yaml`.
First, clone this repository and navigate to the root directory of the repository:

```bash
git clone https://github.com/DendrouLab/panpipes.git
cd panpipes
```

Then, create the conda environment and install the nightly version of panpipes using the following command:

```bash
conda env create --file=pipeline_env.yaml 
conda activate pipeline_env
pip install -e .
```

Oxford BMRC Rescomp users find additional advice on the installation [here](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md).

## Releases

Since `panpipes v.1.2.0` the spatial workflows use SpatialData as the data format for spatial assays. Check out the documentation at https://panpipes-pipelines.readthedocs.io

Since `panpipes v0.4.0`, the `ingest` workflow expects different headers for the RNA and Protein modalities from the previous releases.
Check the example [submission file](https://github.com/DendrouLab/panpipes/blob/main/docs/usage/sample_file_qc_mm.md) and the [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/usage/setup_for_qc_mm.html) for more detailed instructions.

## Citation

[Curion, F., Rich-Griffin, C., Agarwal, D. et al. Panpipes: a pipeline for multiomic single-cell and spatial transcriptomic data analysis. Genome Biol 25, 181 (2024). 
doi: https://doi.org/10.1186/s13059-024-03322-7](https://link.springer.com/article/10.1186/s13059-024-03322-7)


## Contributors

Created by Charlotte Rich-Griffin and Fabiola Curion.
Maintained by Fabiola Curion and Sarah Ouologuem.
Additional contributors: Sarah Ouologuem, Devika Agarwal, Lilly May, Kevin Rue-Albrecht, Giulia Garcia, Wojciech Lason, Lukas Heumos, Jarne Belien.
