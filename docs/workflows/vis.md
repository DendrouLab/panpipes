Visualization
=============

The `vis` workflow can take a `MuData` object from any of the other `Panpipes` workflows and visualize the following (where applicable):

- for custom markers: 
    - heatmaps
    - dotplots
    - embeddings such as PCA and UMAP coloured by the marker expression 
- for categorical variables: 
    - barplots,  e.g the number of cells per sample ID
    - stacked barplots, e.g. the number of cells in each cluster split by sample ID
    - embeddings such as PCA and UMAP coloured by the categorical variable
- for continuous variables: 
    - violin plots split by categorical variables, e.g doublet score per disease group
    - embeddings such as PCA and UMAP coloured by the continuous variable
- for paired markers or metrics: 
    - scatter plots, e.g. scatter plot of `total_counts` on the x axis and `n_genes_by_counts` on the y axis  

For plots with custom markers follow the guidelines in [Gene list formats](../usage/gene_list_format.md) to create the input gene lists.


## Steps to run:
1.  Activate conda environment `conda activate pipeline_env`
2.  Generate the config and log files with `panpipes vis config` and edit the pipeline.yml file
3.  (Optional) Prepare gene list files, [more details here](../usage/gene_list_format.md)
4.  Run complete workflow with `panpipes vis make full --local`

The [Visualization tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/visualization/vis_with_panpipes.html) guides you through the visualization step by step. 


## Expected structure of MuData object
The `vis` workflow can take the MuData outputs from any of the other `Panpipes` workflows, but make sure not to ask the pipeline to produce visualisations using data that has not been computed, e.g. plotting UMAPs using a MuData object which does not contain anything in the X_umap slot.
