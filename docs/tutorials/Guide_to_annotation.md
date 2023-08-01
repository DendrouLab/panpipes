# Guide to annotation for the Dendrou lab pipeline

Tom Thomas and Charlotte Rich Griffin - 13-04-2021

This is essentially a guide to bring together efficient ways of working with the group pipeline output.
It should hopefully also help reduce time spent Googling pandas methods for wrangling the anndata object

## Inputs and adjuncts

inputs: 
- batch corrected n_neighbors h5ad object from the clustering pipeline output directory
- all_res_clusters_list.txt file from the clustering pipeline output directory
       
adjuncts: 

- At the same time as annotating, open up clustree in parallel, and also open up the excel sheet carrying differentially expressed genes
- To get a broad idea of the clusters and how they develop visually on the umap, also open up umap_all_res.png file from the figures folder within the n_neighbour specific clustering pipeline output
- The resolution specific figures folder is also well worth examining as this will enable you to identify straight away the broad clusters (we have included broad gene panels for the various cell fractions) and potential co-variate specific cell abundances

NOTE: In our experience clustree is a simply a guide, interesting biology might occur at varying levels on the clustree output. Annotate across resolutions as biologically indicated

## Aims:
Within the broader Dendrou workflow - the aims of the annotation step is to enable the following:

1. create annotated anndata objects - subsequently combine the annotated fractions to create a fully annotated dataset for the experiment
2. write out annotated barcode (subpopulation, minor, major and bucket), as well as umap co-ordinates for further downstream analyses

BEFORE STARTING: identify broad populations and trends within the dataset by first screening the adjunct pictorial outputs. At this point, edit the clustree with any potential identified broad clusters using the adjunct resources - this will help guide the annotation process. Select three resolutions where there is minimal cross over in the clustree and might be biologically meaningful (this is a guess initially and you can refine as you explore the dataset)


```python
#load necessary packages
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import anndata
import seaborn as sns
```


```python
#set visual settings - this is my own private preference, edit as you wish
```


```python
sc.set_figure_params(scanpy=True, dpi=300, dpi_save=150, frameon=True, vector_friendly=True, fontsize=15, figsize=[35, 35], color_map=None, format='pdf', facecolor=None, transparent=False, ipython_format='png2x')
```


```python
#load the necessary n_neighbours object from the clustering pipeline output
```


```python
adata = sc.read('bcell_drharmony_meteuclidean_neighbors.h5ad')
```


```python
#load cluster data from clustering pipeline output
```


```python
data = pd.read_csv('all_res_clusters_list.txt',sep="\t",index_col=0)
```


```python
#combine both 
```


```python
adata.obs = adata.obs.merge(data, left_index=True, right_index=True)
```


```python
#set the resolutions as identified above of interest as a category
```


```python
adata.obs['leiden_res_0.3'] = adata.obs['leiden_res_0.3'].astype('category')
```

## Identify soupy clusters
You can do this by using the doublet scores (from scrublet) which the pipeline will have added for you.
Mark these on your clustree so that you can forget about these clusters


```python
sns.catplot(x=adata.obs['leiden_res_0.3'], y=adata.obs['doubletscores'], kind = "box", data=adata.obs)
```

## Annotation cell populations
Annotating further from this point on revolves around:

1. reading the top marker genes on the excel sheet and researching them - annotate on the clustree to keep track 
2. comparing expression levels of genes across clusters by visualisation

VISUALISATIONS - I tend to use dotplots as a mainstay for all my annotations, using umap can be helpful to see whether there might be subclusters driven by a single gene within a cluster. Looking for CD3D for ex: helped me distinguish that within our NKT cluster, there was a distinct population which did not express CD3D. On further subclustering, we were able to identify the NK population


```python
#set variables
genes = ['TNF']
```

VISUALISATIONS 1: UMAP + gene on umap (if you have specific queries about a gene)


```python
sc.pl.umap(adata, color=['leiden_res_0.3','TNF'], size = 15, legend_loc = 'on data', legend_fontsize = 'large',add_outline=True)
```

VISUALISATIONS 2: dotplots


```python
genes = ['TNFRSF13B','CXCR5','SELL','CD27','CCR7']
```


```python
sc.pl.dotplot(adata, genes, groupby = 'leiden_res_0.3', dendrogram=True)
```

## Data wrangling. I - subset data according to feature or cell type etc.
Use scenario - this is useful for subsetting cell types for further clustering, and also examining cells with specific covariates (inflamed vs non-inflamed etc.) further


```python
#Ensure cluster column is category
```


```python
adata.obs['leiden_res_0.8_revised'] = adata.obs['leiden_res_0.8_revised'].astype('category')
```

subset according to column, and value in column - this isolates cells with value 'Inflamed' from the anndata object using the adata.obs column 'Inflammation'
you can use the same concept to isolate cell types of interest from the annotation column


```python
infl = adata[adata.obs['Inflammation'] == 'Inflamed']
```

## Data wrangling. II - retrieving annotations at multiple levels on the clustree 
Use scenario: Suppose you think that the cluster resolution at 0.8 satisfies the vast majority of your annotation, but you discern that there is a biologically relevant subtype at resolution 1.1
Use this to bring that cluster at res 1.1 to res 0.8
This is really just an exercise in pandas dataframe wrangling


```python
##Ensure both are data type category first
```


```python
adata.obs['leiden_res_0.8'] = adata.obs['leiden_res_0.8'].astype('category')
adata.obs['leiden_res_1.1'] = adata.obs['leiden_res_1.1'].astype('category')
```


```python
##Isolate both columns into a single dataframe
```


```python
df_test = pd.DataFrame(data={'Cell': list(adata.obs.index), 'leiden_res_1.1': list(adata.obs['leiden_res_1.1']), 'leiden_res_0.8': list(adata.obs['leiden_res_0.8'])})
df_test.set_index('Cell',inplace=True)
```


```python
##Anywhere, that resolution 1.1 deems to be 12 (your cluster of interest), replace at leiden 0.8 with 19 (new cluster that you have deemed necessary on 0.8 res)
```


```python
df_test.loc[df_test['leiden_res_1.1'] == 12, "leiden_res_0.8"] = 19
```


```python
##double check that this has worked as you intended
```


```python
df_test["leiden_res_0.8"].unique()
```


```python
##delete now defunct resolution 1.1
```


```python
del df_test['leiden_res_1.1']
```


```python
##rename column to keep tracking easier
```


```python
df_test.columns = ['leiden_res_0.8_revised']
```


```python
##merge this column into adata.obs section
```


```python
adata.obs = adata.obs.merge(df_test, left_index=True, right_index=True)
```


```python
##explictly put this as a category and then treble check that you wanted to happen, actually happened on the umap
```


```python
adata.obs['leiden_res_0.8_revised'] = adata.obs['leiden_res_0.8_revised'].astype('category')
```


```python
sc.pl.umap(adata, color=['leiden_res_0.8_revised'], size = 15, legend_loc = 'on data', legend_fontsize = 'large',add_outline=True)
```

## Data wrangling.III Merge clusters at resolution of your choice

Ensure both are data type category first

Use scenario: Suppose you think that at the resolution of your choice, there is no sufficient biological distinction between two clusters, and you would rather merge them
Use this to bring that cluster at res 1.1 to res 0.8
This is really just an exercise in pandas dataframe wrangling

As before in Data wrangling - II, I like to create a separate dataframe before merging it into the anndata object, once I am happy with annotations.
This is a personal choice - it just means I am not messing up and having to restart if there are any issues in the code.


```python
df_test = pd.DataFrame(data={'Cell': list(adata.obs.index), 'leiden_res_0.9': list(adata.obs['leiden_res_0.9'])})
df_test.set_index('Cell',inplace=True)
```

Checking the per cluster number of cells is good, as it allows you to sanity check whether this operation has been conducted as expected


```python
df_test['leiden_res_0.9'].value_counts()
```

Merging cluster 10 and 0: Any cell that is annotated as 10, make cluster 0


```python
df_test.loc[df_test['leiden_res_0.9'] == 10, "leiden_res_0.9"] = 0
```


```python
#Sanity check to see if you see any more cells annotated as cluster 10 #there shouldn't be!
```


```python
df_test["leiden_res_0.9"].unique()
```


```python
#Sanity check numbers count to ensure those cells have been added to cluster 0
```


```python
df_test['leiden_res_0.9'].value_counts()
```


```python
#rename column - as otherwise when you merge, two columns will have the same name so pandas will append _x and _y to both column names
```


```python
df_test.columns = ['leiden_res_0.9_revised']
```


```python
#carry out merge
```


```python
adata.obs = adata.obs.merge(df_test, left_index=True, right_index=True)
```


```python
#explictly put this as a category and then treble check that you wanted to happen, actually happened on the umap
```


```python
adata.obs['leiden_res_0.9_revised'] = adata.obs['leiden_res_0.8_revised'].astype('category')
```


```python
sc.pl.umap(adata, color=['leiden_res_0.9_revised'], size = 15, legend_loc = 'on data', legend_fontsize = 'large',add_outline=True)
```

## Data wrangling.IV Rename annotations


```python
#check that the resolution of interest is a category data type
```


```python
adata.obs['leiden_res_0.3'].cat.categories
```


```python
#remove soupy cells from final anndata object
```


```python
adata2 = adata[adata.obs['leiden_res_0.3'] != 5]
adata2 = adata2[adata2.obs['leiden_res_0.3'] != 6]
adata2 = adata2[adata2.obs['leiden_res_0.3'] != 10]
adata2 = adata2[adata2.obs['leiden_res_0.3'] != 11]
adata2 = adata2[adata2.obs['leiden_res_0.3'] != 8]
adata2 = adata2[adata2.obs['leiden_res_0.3'] != 3]
```


```python
#double check that the objects have been removed
```


```python
adata2.obs['leiden_res_0.3'].cat.categories
```


```python
#rename annotations
```


```python
adata2.rename_categories('leiden_res_0.3', ['Memory B cell','Follicular B cell','TNF hi memory cell','FCRL4 mem B cell','IFN hi mem B cell','IgG+ mem B cell'])
```


```python
#delete other clusters
```


```python
del adata2.obs['leiden_res_0.1']
del adata2.obs['leiden_res_0.2']
del adata2.obs['leiden_res_0.4']
del adata2.obs['leiden_res_0.5']
del adata2.obs['leiden_res_0.6']
del adata2.obs['leiden_res_0.7']
del adata2.obs['leiden_res_0.8']
del adata2.obs['leiden_res_0.9']
del adata2.obs['leiden_res_1']
del adata2.obs['leiden_res_1.1']
del adata2.obs['leiden_res_1.2']
del adata2.obs['leiden_res_1.3']
del adata2.obs['leiden_res_1.4']
```


```python
#rename the annotation column
```


```python
adata2.obs.rename(columns={"leiden_res_0.3": "final"})
```

## Final marker visualisation 
gene expression, percentage of cells expressing and cell frequency
make list of genes that could be put into a paper


```python
genes = ['TNFRSF13B','CXCR5','SELL','CD27','CCR7','IGHD','FCER2','TCL1A', 'TNF','FCRL4']
```


```python
a = sc.pl.dotplot(adata2, genes, groupby = 'final', return_fig=True)
```


```python
a.add_totals(color='#0072B2', sort=True).style(dot_edge_color='black', dot_edge_lw=0.25).show()
```


```python
#write out final anndata object
```


```python
adata2.write('bcell_annotated.h5ad')
```

once combined - write out all annotations plus the umap co-ordinates


```python
#annotations
```


```python
df_distribution = pd.DataFrame(data={'barcode': list(adata2.obs.index), 'cell_type': list(adata2.obs['final'])})
df_distribution.set_index('barcode',inplace=True)
```


```python
#umap
```


```python
a = pd.DataFrame(adata.obs.index, columns = ["barcode"])
```


```python
b = pd.DataFrame(adata.obsm['X_umap'], columns = ['X1','X2'])
```


```python
res = a.join(b)
res = res.set_index('barcode')
```


```python
#export annotations and umap
```


```python
export = df_distribution.obs.merge(res, left_index=True, right_index=True)
```


```python
export.to_csv('export.tsv', sep='\t')
```
