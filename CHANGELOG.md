
# Releases

## [latest]

### added 
- moved from MuData to SpatialData 
- xenium ingestion & ingest_xenium GitHub action
- `export_gene_by_spot` for Cell2Location
- separate sample submission files for the different spatial technologies 
- separate GitHub actions for Cell2Location & Tangram 
- separate GitHub actions for MERSCOPE & MERFISH


### fixed


### dependencies
- pinned "spatialdata==0.2.6", "spatialdata-io==0.1.6", "dask==2024.12.1" as temporary fix


## v1.1.0

### added 

### fixed

- fixed calls to new imports matplotlib 
- fixed scanpy/muon latest
- changed default adversarial_training = False to ensure multivi runs with scvi-tools 1.1.3 https://github.com/scverse/scvi-tools/issues/2581
- updated integration01 action to use updated scvi-tools version

### dependencies

- matplotlib 3.9.0 


## v1.0.0

### added
- added documentation on ingesting custom h5ad/h5mu objects

### fixed


### dependencies

- pinned matplotlib 3.8 as temporary fix 
 next release will unpin and fix get_cmap ([API](https://matplotlib.org/stable/api/cm_api.html#:~:text=The%20colormap%20name%20can%20then%20be%20used%20as%20a%20string%20argument%20to%20any%20cmap%20parameter%20in%20Matplotlib.%20It%20is%20also%20available%20in%20pyplot.get_cmap.)) and legendHandles ([API](https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.9.0.html))


## v0.5.0

### added

- added Tangram to `deconvolution_spatial`
- added scib metrics calculation to `integration`, using the [scib-metrics package](https://scib-metrics.readthedocs.io/en/latest/index.html)
- [extra documentation](https://panpipes-pipelines.readthedocs.io/en/latest/)

### fixed

- fixed error in `vis`
  - error occurred when only wanting to plot continuous or categorical variables (or neither), not both
- fixed error in `refmap`
  - high threads was not recognised, now fixed.

### dependencies
- All the dependencies have been updated. 
  - Python>=3.10 required
- added seeds to all scvi tasks for reproducibility

## v0.4.1

### added

- added example multiome submission file 10X_h5
- added example multiome submission file cellranger
- workflows & tutorials for `qc_spatial`, `preprocess_spatial`, and `deconvolution_spatial` to readthedocs
- tutorial for `vis`
- added PCA parameters in pipeline_preprocess.py for PROT modality to fix issue #120
- added full control of dimred params for all modalities in pipeline_preprocess.py
- more info on custom genes format files added to documentation
- parsing summary files for cellranger multi version < 7
- added checks for n_pcs in run_neighbors_method_choice
- added filtering by HVF for atac

### fixed

- changed typo in tutorial paths for clustering and deconvolution
- fix io to read cellranger outs folder for atac.
- fixes to refmap workflow
- typos & capitalization in the pipeline.yml files of `qc_spatial`, `preprocess_spatial`, and `deconvolution_spatial`, `vis`
- remove `assay`, `sample_prefix`, and `modalities` parameters from the `qc_spatial` pipeline.yml
- remove `sample_prefix` and `modalities` parameters from the `preprocess_spatial` pipeline.yml
- fixed error in `preprocess_spatial` when `filtering: run: False`
    -> now able to run no filtering without needing to save the MuData in `filtered.data` before running the pipeline
- fixed error in `vis`
  - change PARAMS['custom_markers_minimal'] -> PARAMS['custom_markers']['files']['minimal']
- fix to avoid rerunning HVF and explicitly check X layer before normalization in pipeline_preprocess.py
- fix plotting of umaps after batch correction
- fix fetching string scvi if present in mudata for wnn
- fixed lsi requirement for atac
- fixed top features for atac
- fixed filtering HVG for rna
- moved pynndescent to PyPi dependencies

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

## v0.2.0

- First public version of panpipes
- contains qc_mm, preprocess, intergration, clustering
