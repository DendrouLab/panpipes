Running pipeline modules from different entry points.
=====================================================

For circumstances where you already have a qc'd Mudata object, you can run any of the other pipelines. In order to do so, you need to follow some specific file naming conventions.
You must always have a column called sample_id within each modality in the mudata object, as well as in the top `mdata.obs`. 

## Preprocess

It is possible to run the Preprocess pipeline starting from one combined mudata object. It must contain raw data in each modalities X slots.  

If your data is pre-filtered call your anndata object `[PROJECT_PREFIX]_filt.h5mu` and set `filtering_run: False` in the `pipeline.yml`.

If you have not filtered your data then you can set run filtering_run: True, and set the remaining parameters, but make sure you review how the dynamic filtering works [here](filter_dict_instructions).

## Integration 

Coming soon.

<!-- Raw counts in the layers, normalised/scaled counts in the X -->

## Clustering

To run clustering_scanpy without the prior steps, you will need to produce a
[PROJECT_PREFIX].h5mu 
It should contain: 
- log normalised data in the `mdata['rna'].X` slot
- highly variable genes calculated
- normalised data in the other modality X slots.
- pca caluclated for each modality (using scanpy's `sc.pp.pca` function)


## Refmap

Coming soon.
<!-- Raw counts in the layers -->


## Vis

Any h5mu object but containing the columns and normalised layers that you expect to plot.
