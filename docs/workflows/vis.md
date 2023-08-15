Visualisation
=============

`panpipes visualisation` can take a MuData object from any of the other `panpipes` pipelines and visualise the following (where applicable).
- custom markers on heatmaps/dotplots
- custom markers on embeddings such as PCA and UMAP
- barplots on categorical variables . e.g the number of cells per sample_id
- stacked barplots on categorical variables e.g. the number of cells in each cluster split by sample_id
- violin plots of continuous variables split by categorical variables e.g doublet score per disease group

For plots with custom markers follow the guidelines in [Gene list formats](../usage/gene_list_format.md) to create the input gene lists.


## Steps:
1.  In a new folder, generate config file for integration,
    `panpipes vis config` and edit the pipeline.yml file.
2.  Prepare plotting gene list files
    -   [more
        details](https://github.com/DendrouLab/panpipes/blob/main/docs/gene_list_format.md)
3.  Run complete refmap pipeline with `panpipes vis make full`

To repeat the pipeline after editing the pipeline.yml, delete the files
in the log folder and repeat step 3.

## Expected structure of MuData object
`Panpipes vis` can take the MuData outputs from any of the other `panpipes` workflows, but make sure not to ask the pipeline to produce visualisations using data that has not been computed, e.g. plotting UMAPs using a MuData object which does not contain anything in the X_umap slot.
