# Panpipes - multimodal single cell pipelines 

Created and Maintained by Charlotte Rich-Griffin and Fabiola Curion  
Additional contributors: Sarah Ouologuem, Devika Agarwal and Tom Thomas 

See out [documentation](https://panpipes-pipelines.readthedocs.io/en/latest/) and our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.11.532085v1):  
Panpipes: a pipeline for multiomic single-cell data analysis  
Charlotte Rich-Griffin*, Fabiola Curion*, Tom Thomas, Devika Agarwal, Fabian J. Theis, Calliope A. Dendrou.  
bioRxiv 2023.03.11.532085;  
doi: https://doi.org/10.1101/2023.03.11.532085


# Introduction
These pipelines use cgat-core pipeline software

Available pipelines:
1. "ingest" : for the ingestion of data and computation of QC metrics' 
2. "preprocess" : for filtering and normalising of each modality
3. "integration" : integrate and batch correction using  single and multimodal methods
4. "clustering" : cell clustering on single modalities
5. "refmap" : transfer scvi-tools models from published data to your data
6. "vis" : visualise metrics from other pipelines in context of experiment metadata
7. "qc_spatial": for the ingestion of Spatial Transcriptomics data (vizgen, visium) and computation of QC metrics
8. "preprocess_spatial" : for filtering and normalising of ST data
9. "deconvolution": for the cell type deconvolution of ST slides


# Installation and configuration
See [installation instructions here](https://panpipes-pipelines.readthedocs.io/en/latest/install.html)


Oxford BMRC Rescomp users find additional advice in [docs/installation_rescomp](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md)

