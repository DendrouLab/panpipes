

## [main]
### added

### fixed

### dependencies

## v0.3.1
- set default matplotlib<=3.7.3 to avoid issue #104. 
- add workflow for `qc_spatial` to readthedocs
- remove `assay` and `sample_prefix` parameters from the `qc_spatial` pipeline.yml 
- add workflow for `preprocess_spatial` to readthedocs

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

