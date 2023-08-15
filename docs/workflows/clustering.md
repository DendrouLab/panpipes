Clustering
==========

The panpipes clustering pipeline runs the following steps:
- detect or generate nearest neighbor graph - In order to take advantage of batch corrected data, we recommend using the `integration` pipeline which generates nearest neighbor graphs based on batch corrected/integrate data. If you wish to create a new nearest neighbor graph (perhaps to test a different parameter set), then fill in the yml as follows, and repeat for each modality: 
```
neighbors:
  rna:
    use_existing: False
    # number of Principal Components to calculate for neighbours and umap:
    dim_red: X_harmony (or any reduced dimensional represenation in your MuData object)
    #how many components to use for clustering
    n_dim_red: 30
    # number of neighbours
    k: 30
    # metric: euclidean | cosine
    metric: euclidean
    # scanpy | hnsw (from scvelo)
    method: scanpy
```
- generate UMAPS with different minimum distance parameters
- calculate leiden or louvain clustering, at a range of resolutions, saving them as csv and in the MuData object.
- visualise clustering on UMAPs and as a [clustree diagram](https://lazappi.github.io/clustree/articles/clustree.html)
- identify top marker genes for each cluster at each clustering resolution, using [`scanpy.tl.rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) saved as csvs. 


## Steps to run:
1.  In a new folder, generate config file for integration,
    `panpipes clustering config` and edit the pipeline.yml file.
2.  Run the clustering pipeline
    `panpipes clustering make cluster_analysis`. This will do the
    initial nearest neighbours and clustering for the parameters you
    specify.
3.  Decide on the best values for k nearest neighbours based on UMAPs
    and clustree results. (Optionally) Once decided delete the folders for the
    parameters you don't need and delete those from the pipeline.yml.
4.  Find markers for each of your cluster resolutions with
    `panpipes clustering make marker_analysis` (Again you could run all
    the clustering pipeline at once with `panpipes clustering make full`
    but by making decisions along the way you'll reduce the computation
    and file size burden of the pipeline)

## Expected structure of MuData object
The ideal way to run `panpipes clustering` is to use the output mudata file from `panpipes preprocess` or `panpipes integration`, as this will make sure the MuData object has correctly names layers and slots. 

The bare minimum MuData object required is normalised data in the X slot of each modality,  a sample_id column in each slot of the obs and the outer obs, and the layers indicated in the find_markers section of the yml.