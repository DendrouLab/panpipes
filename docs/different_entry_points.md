# Running pipeline modules from different entry points.

For circumstances where you already have a qc'd Mudata object
You must always have a column called sample_id which groups the data in some way. 

## Preprocess
It is possible to run the Preprocess pipeline starting from one combined anndata object containing all your samples containing raw counts in the X slot, 
either with or without running filtering first.

If your data is pre-filtered call your anndata object [PROJECT_PREFIX]_filt.h5mu and set filtering_run: False.

If you have not filtered your data then you can set run filtering_run: True, and set the remaining parameters, but make sure you review how the dynamic 

## Integration 
Coming soon.


## Clustering
To run clustering_scanpy without the prior steps, you will need to produce a
[PROJECT_PREFIX].h5mu 
It should contain: 
- log normalised data in the mdata['rna'].X slot
- highly variable genes calculated
- normalised data in the other modality X slots.
- pca 


## Refmap

Coming soon.


## Vis
- Any h5mu object but containing the columns and noramlised layers that you expect to plot.
