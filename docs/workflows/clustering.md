# Clustering

The panpipes clustering pipeline runs the following steps:

- Detect or generate nearest neighbor graph 
- Generate UMAPS with different minimum distance parameters
- Calculate leiden or louvain clustering, at a range of resolutions, saving them as csv and in the MuData object.
- Visualise clustering on UMAPs and as a [clustree diagram](https://lazappi.github.io/clustree/articles/clustree.html). Note that it is important to specify more than 1 resolution for clustree to work.
- Identify top marker genes for each cluster at each clustering resolution, using [`scanpy.tl.rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) saved as csvs.



## Steps to run



1. In a new folder, generate config file for clustering,
    `panpipes clustering config` and edit the pipeline.yml file.
2. Run the clustering pipeline
   you can run the full pipeline with `panpipes clustering make full` or follow this step by step guide:

1.  `panpipes clustering make cluster_analysis`. This will do the
    initial nearest neighbours and clustering for the parameters you
    specify.
2. Decide on the best values for k nearest neighbours based on UMAPs
    and clustree results. (Optionally) Once decided delete the folders for the
    parameters you don't need and delete those from the pipeline.yml.
3. Find markers for each of your cluster resolutions with
    `panpipes clustering make marker_analysis` (Again you could run all
    the clustering pipeline at once with `panpipes clustering make full`
    but by making decisions along the way you'll reduce the computation
    and file size burden of the pipeline)

## Expected structure of MuData object

The ideal way to run `panpipes clustering` is to use the output MuData file from `panpipes preprocess` or `panpipes integration`, as this will make sure the MuData object has correctly names layers and slots.

The bare minimum MuData object required is normalised data in the `.X` slot of each modality, a sample_id column in each slot of the obs and the outer obs, and the layers indicated in the find_markers section of the yml.

Check the [clustering tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/clustering/clustering_tutorial.html) for an example of how to run the workflow.


## A note on neighbors graph:
In order to take advantage of batch corrected data, we recommend using the `integration` workflow which generates K nearest neighbor graphs (KNN) based on batch corrected/integrated data, before running `clustering`. If you wish to create a new nearest neighbor graph (perhaps to test a different parameter set), then make sure to set the `use_existing: False` in the clustering `yml` configuration file and fill in the yml with your custom params, and repeat for each modality for which you want to recompute the neighbors :

```yaml
neighbors:
  rna:
    use_existing: False
    # which representation in .obsm to use for nearest neighbors
    # if dim_red=X_pca and X_pca not in .obsm, will be computed with default parameters
    dim_red: X_pca
    #how many components to use for clustering
    n_dim_red: 30
    # number of neighbours
    k: 30
    # metric: euclidean | cosine
    metric: euclidean
    # scanpy | hnsw (from scvelo)
    method: scanpy
```


