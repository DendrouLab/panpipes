
# Spatial Deconvolution YAML

In this documentation, the parameters of the `deconvolution_spatial` yaml file are explained. 
This file is generated running `panpipes deconvolution config`.
In general, the user can leave parameters empty to use defaults. <br> The individual steps run by the pipeline are described in the [spatial deconvolution workflow](../workflows/deconvolute_spatial.md). 

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)


## 0. Compute Resource Options

| `resources` |  |
| --- | --- |
| `threads_high` | __`int`__ (default: 1) <br> Number of threads used for high intensity computing tasks. |
| `threads_medium` | __`int`__ (default: 1) <br> Number of threads used for medium intensity computing tasks. For each thread, there must be enough memory to load your mudata and do computationally light tasks. |
| `threads_low` | __`int`__ (default: 1) <br> Number of threads used for low intensity computing tasks. For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.|

|  |  |
| ---- | --- |
| `condaenv` | __`str`__ (default: None) <br> Path to conda environment that should be used to run `Panpipes`. Leave blank if running native or your cluster automatically inherits the login node environment. |


## 1. Input Options
With the `deconvolution_spatial` workflow, one or multiple spatial slides can be deconvoluted in one run. For that, a `MuData` object for each slide is expected, with the spatial data saved in `mdata.mod["spatial"]`. The spatial slides are deconvoluted **using the same reference**. For the reference, one `MuData` with the gene expression data saved in `mdata.mod["rna"]` is expected as input. Please note, that the same parameter setting is used for each slide. <br> For the **spatial** input, the workflow, therefore, reads in **all `.h5mu` objects of a directory** (see below). **The spatial and single-cell data thus need to be saved in different folders.**

| `input` |  |
| ---- | --- |
| `spatial` | __`str`__ (not optional) <br> Path to folder containing one or multiple `MuDatas` of spatial data. The pipeline is reading in all `MuData` files in that folder and assuming that they are `MuDatas` of spatial slides.|
| `singlecell` | __`str`__ (not optional) <br> Path to the MuData **file** (not folder) of the reference single-cell data.|

## 2. Cell2Location Options

For each deconvolution method you can specify whether to run it or not: 
| |  |
| ---- | --- |
| `run` | __`bool`__ (default: None) <br> Whether to run Cell2location|

### Feature Selection 

You can select genes that are used for deconvolution in two ways. The first option is to provide a reduced feature set as a csv-file that is then used for deconvolution. The second option is to perform gene selection [according to Cell2Location](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html). <br> Please note, that gene selection is **not optional**. If no csv-file is provided, feature selection [according to Cell2Location.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html) is performed. 


| `feature_selection` |  |
| ---- | --- |
| `gene_list` | __`str`__ (default: None) <br>  Path to a csv file containing a reduced feature set. A header in the csv is expected in the first row. All genes of that gene list need to be present in both, spatial slides and scRNA-Seq reference.|
| `remove_mt` | __`bool`__ (default: True) <br> Whether to remove mitochondrial genes from the dataset. This step is performed **before** running gene selection. |
| `cell_count_cutoff` | __`int`__ (default: 15) <br> All genes detected in less than cell_count_cutoff cells will be excluded. Parameter of the [Cell2Location's gene selection function.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html)|
| `cell_percentage_cutoff2` | __`float`__ (default: 0.05) <br> All genes detected in at least this percentage of cells will be included. Parameter of the [Cell2Location's gene selection function.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html)|
| `nonz_mean_cutoff` | __`float`__ (default: 1.12) <br> Genes detected in the number of cells between the above-mentioned cutoffs are selected only when their average expression in non-zero cells is above this cutoff.  Parameter of the [Cell2Location's gene selection function.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html) | 


### Reference Model

| `reference` |  |
| ---- | --- |
| `labels_key` | __`str`__ (default: None) <br> Key in `.obs` for label (cell type) information. |
| `batch_key` | __`str`__ (default: None) <br> Key in `.obs` for batch information. |
| `layer` | __`float`__ (default: None) <br> Layer in `.layers` to use for the reference model. If None, `.X` will be used. Please note, that Cell2Location expects raw counts as input.|
| `categorical_covariate_keys` | __`str`__ (default: None) <br> Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space).| 
| `continuous_covariate_keys` | __`str`__ (default: None) <br> Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space)| 
| `max_epochs` | __`int`__ (default: _np.min([round((20000 / n_cells) * 400), 400])_) <br> Number of epochs.| 
| `use_gpu` | __`bool`__ (default: True) <br> Whether to use GPU for training. | 


### Spatial Model 


| `spatial` |  |
| ---- | --- |
| `batch_key` | __`str`__ (default: None) <br> Key in `.obs` for batch information. |
| `layer` | __`float`__ (default: None) <br> Layer in `.layers` to use for the reference model. If None, `.X` will be used. Please note, that Cell2Location expects raw counts as input.|
| `categorical_covariate_keys` | __`str`__ (default: None) <br> Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space).| 
| `continuous_covariate_keys` | __`str`__ (default: None) <br> Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space)| 
| `N_cells_per_location` | __`int`__ (not optional) <br> Expected cell abundance per voxel. Please refer to the [Cell2Location documentation](https://cell2location.readthedocs.io/en/latest/index.html) for more information. | 
| `detection_alpha` | __`float`__ (not optional) <br> Regularization of with-in experiment variation in RNA detection sensitivity. Please refer to the [Cell2Location documentation](https://cell2location.readthedocs.io/en/latest/index.html) for more information.  | 
| `max_epochs` | __`int`__ (default: _np.min([round((20000 / n_cells) * 400), 400])_) <br> Number of epochs.| 
| `use_gpu` | __`bool`__ (default: True) <br> Whether to use GPU for training. | 


###
<br>
You can specify whether both models should be saved with the following parameter: 

| |  |
| ---- | --- |
| `save_models` | __`bool`__ (default: False) <br> Whether to save the reference & spatial mapping models|




## 3. Tangram Options

For each deconvolution method you can specify whether to run it or not: 
| |  |
| ---- | --- |
| `run` | __`bool`__ (default: None) <br> Whether to run Tangram|


### Feature Selection

You can select genes that are used for deconvolution in two ways. The first option is to provide a reduced feature set as a csv-file that is then used for deconvolution. The second option is to perform gene selection via [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) **on the reference scRNA-Seq data**, as [suggested by Tangram](https://tangram-sc.readthedocs.io/en/latest/tutorial_sq_link.html#Pre-processing). The top `n_genes` of each group make up the reduced gene set. <br> Please note, that gene selection is **not optional**. If no csv-file is provided, feature selection via [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) is performed. 


| `feature_selection` |  |
| ---- | --- |
| `gene_list` | __`str`__ (default: None) <br>  Path to a csv file containing a reduced feature set. A header in the csv is expected in the first row. All genes of that gene list need to be present in both, spatial slides and scRNA-Seq reference.|

___Parameters for `scanpy.tl.rank_genes_groups` gene selection___ 
| `rank_genes` |  |
| ---- | --- |
| `labels_key` | __`str`__ (default: None) <br> Which column in `.obs` of the reference to use for the `groupby` parameter of [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) .|
| `layer` | __`str`__ (default: None) <br> Which layer of the reference to use for [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html). If None, `.X` is used.|
| `n_genes` | __`int`__ (default: 100) <br> How many top genes to select of each `groupby` group|
| `test_method` | __`str ['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']`__ (default: 't-test_overestim_var') <br> Which test method to use.|
| `correction_method` | __`str ['benjamini-hochberg', 'bonferroni']`__ (default: ' benjamini-hochberg') <br> Which p-value correction method to use. Used only for 't-test', 't-test_overestim_var', and 'wilcoxon'. |

### Model

| `model` |  |
| ---- | --- |
| `labels_key` | __`str`__ (default: None) <br> Key in `.obs` for label (cell type) information. |
| `num_epochs` | __`int`__ (default: 1000) <br> Number of epochs. | 
| `device` | __`str`__ (default: 'cpu') <br> Which device to use. | 
| `kwargs` | In `kwargs`, the user has the possibility to specify parameters for [tangram.mapping_utils.map_cells_to_space](https://tangram-sc.readthedocs.io/en/latest/classes/tangram.mapping_utils.map_cells_to_space.html?highlight=mapping_utils%20map_cells_to_space#tangram.mapping_utils.map_cells_to_space). You can add or remove any parameters of the function.| 

