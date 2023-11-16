

## [main]

### added
- added example multiome submission file 10X_h5
- added example multiome submission file cellranger
- workflows & tutorials for `qc_spatial`, `preprocess_spatial`, and `deconvolution_spatial` to readthedocs

### fixed

- changed typo in tutorial paths for clustering and deconvolution
- fix io to read cellranger outs folder for atac. 
- typos & capitalization in the pipeline.yml files of `qc_spatial`, `preprocess_spatial`, and `deconvolution_spatial`
- remove `assay`, `sample_prefix`, and `modalities` parameters from the `qc_spatial` pipeline.yml 
- remove `sample_prefix` and `modalities` parameters from the `preprocess_spatial` pipeline.yml
- fixed error in `preprocess_spatial` when `filtering: run: False`
    -> now able to run no filtering without needing to save the MuData in `filtered.data` before running the pipeline

### dependencies

## v0.4.0
Big Change! the submission files for the `ingest` workflow have now changed! we require the paths to the Gene expression (RNA/GEX) and Protein (ADT) to have the following headers.


| sample_id | rna_path    | rna_filetype | prot_path    | prot_filetype |
| --------- | ----------- | ------------ | ------------ | ------------- |
| sampleX   | path/to/rna | 10X_h5       | path/to/prot | 10x_h5        |
|           |             |              |              |               |


See tutorials for examples of submission files.



### added
- merged PR #111:
  - LSI in panpipes_preprocess is run on the highly variable features
  - n_comp for LSI
### fixed
- changed all instances of ADT into PROT
- changed all instances of GEX to RNA
- changed the params to fix plotting as mentioned in issue #41
- typo in readme
- set default seaborn <=0.12.2 to avoid issue #104, #126

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


