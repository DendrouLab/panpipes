# Running pipeline modules from different entry points.

For circumstances where you already have a quality controlled Mudata object, you can run any of the other pipelines. In order to do so, you need to follow some specific file naming conventions.
You must always have a column called `sample_id` within the `.obs` of each modality in the MuData object, as well as in the top `mdata.obs`. Please note that the `anndata`/`MuData` object must have the gene/feature name as the index of the `.var` dataframe to be ingested in panpipes. 

## Preprocess

It is possible to run the Preprocess pipeline starting from one combined MuData object. It must contain raw data in each modalities X slots.  

If your data is pre-filtered call your anndata object `[PROJECT_PREFIX]_filt.h5mu` and set `filtering_run: False` in the `pipeline.yml`.

If you have not filtered your data then you can set run filtering_run: True, and set the remaining parameters, but make sure you review how the dynamic filtering works [here](filter_dict_instructions).

For information on normalization methods included in panpipes please check the [Normalization methods](./Normalization%20methods.md) section.

## Integration

You can supply a MuData object with normalized data in the `.X` slot of every modality directly to the `integration` workflow.

If you have previously computed latent representations for the modalities in the MuData, these should be in the relevant `.obsm` slot.
If no dimensionality reduction is present, the `integration` workflow will default PCA for every modality, with [default parameters](https://github.com/scverse/scanpy/blob/master/scanpy/tools/_utils.py#L28) adjusting for number of pcs to be `adata.n_vars -1`

## Clustering

To run clustering_scanpy without the prior steps, you will need to produce a
[PROJECT_PREFIX].h5mu
for the modalities you want to analyse, it should contain:

- normalized data (specify the normalized layer name in the pipeline.yml)

Desireable:

- highly variable features
- dimensionality reduction

If no dimensionality reduction and HVF are present, the pipeline defaults to calculate PCA with default parameters on the full set of available features, and will use the knn graph built on the PCA dimensions to run clustering.

## Refmap

Refmap accepts as query any mudata with raw counts in the `mdata['rna'].X` layer, the reference model and, optionally, the reference AnnData object.
For more info on formatting the query see the [refmap workflow](https://panpipes-pipelines.readthedocs.io/en/latest/workflows/refmap.html)

## Vis

the Vizualization (Viz) workflow accepts any .h5mu object containing the columns and normalised layers that you expect to plot.
