

## [main]
### added
- workflows & tutorials for `qc_spatial`, `preprocess_spatial`, and `deconvolution_spatial` to readthedocs

### fixed
- remove `assay` and `sample_prefix` parameters from the `qc_spatial` pipeline.yml 
- remove `sample_prefix` and `modalities` parameters from the `preprocess_spatial` pipeline.yml
- fixed error in `preprocess_spatial` when `filtering: run: False`
    -> now able to run no filtering without needing to save the MuData in `filtered.data` before running the pipeline

### dependencies

## v0.3.1
- set default matplotlib<=3.7.3 to avoid issue #104. 

## v0.3.0
### added
- Spatial data analysis is now included in panpipes
    - panpipes qc_spatial
    - panpipes preprocess_spatial
    - panpipes_deconvolution_spatial

### fixed
- make sure columns from individual modalities that are not in the multimodal outer obs can be used to

### dependencies
- additional dependencies: squidpy, cell2location, openpyxl

# v 0.2.0
- First public version of panpipes
    - contains qc_mm, preprocess, intergration, clustering

