[![PyPI version](https://badge.fury.io/py/panpipes.svg)](https://badge.fury.io/py/panpipes)
# Panpipes - multimodal single cell pipelines 

Created and Maintained by Charlotte Rich-Griffin and Fabiola Curion  
Additional contributors: Sarah Ouologuem, Devika Agarwal, and Tom Thomas 


**See our [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/) and our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.11.532085v1)**:  
Panpipes: a pipeline for multiomic single-cell data analysis  
Charlotte Rich-Griffin*, Fabiola Curion*, Tom Thomas, Devika Agarwal, Fabian J. Theis, Calliope A. Dendrou.  
bioRxiv 2023.03.11.532085;  
doi: https://doi.org/10.1101/2023.03.11.532085


# Introduction
These workflows use cgat-core pipeline software

Available workflows:
1. "ingest" : for the ingestion of data and computation of QC metrics 
2. "preprocess" : for filtering and normalising of each modality
3. "integration" : integrate and batch correction using  single and multimodal methods
4. "clustering" : cell clustering on single modalities
5. "refmap" : transfer scvi-tools models from published data to your data
6. "vis" : visualize metrics from other pipelines in the context of experiment metadata
7. "qc_spatial" : for the ingestion of spatial transcriptomics (ST) data (Vizgen, Visium) and computation of QC metrics
8. "preprocess_spatial" : for filtering and normalizing ST data
9. "deconvolution_spatial" : for the cell type deconvolution of ST slides


# Installation and configuration
See [installation instructions here](https://panpipes-pipelines.readthedocs.io/en/latest/install.html)


Oxford BMRC Rescomp users find additional advice in [docs/installation_rescomp](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md)

# Releases

`panpipes v0.4.0` is out [now](./CHANGELOG.md)! 

The `ingest` workflow now expects different headers for the RNA and Protein modalities.
Check the example [submission file](https://github.com/DendrouLab/panpipes/blob/main/docs/usage/sample_file_qc_mm.md) and the [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/usage/setup_for_qc_mm.html) for more detailed instructions.

