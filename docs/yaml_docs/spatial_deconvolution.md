<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# Spatial Deconvolution YAML

In this documentation, the parameters of the `deconvolution_spatial` configuration yaml file are explained. 
This file is generated running `panpipes deconvolution_spatial config`.  <br> The individual steps run by the pipeline are described in the [spatial deconvolution workflow](../workflows/deconvolute_spatial.md). 

When running the deconvolution workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.
However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).
You can download the different deconvolution pipeline.yml files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes deconvolution_spatial config`: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_deconvolution_spatial/pipeline.yml)
- `pipeline.yml` file for [Deconvoluting spatial data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/deconvolution/deconvoluting_spatial_data_with_panpipes.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/deconvolution/pipeline.yml)

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)


## 0. Compute Resource Options

<span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.
Specified by the following three parameters:
  - <span class="parameter">threads_high</span> `Integer`, Default: 1<br>
        Number of threads used for high intensity computing tasks. 

  - <span class="parameter">threads_medium</span> `Integer`, Default: 1<br>
        Number of threads used for medium intensity computing tasks.
        For each thread, there must be enough memory to load your SpatialData and do computationally light tasks.

  - <span class="parameter">threads_low</span> `Integer`, Default: 1<br>
  	    Number of threads used for low intensity computing tasks.
        For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.

<span class="parameter">condaenv</span> `String`<br>
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment



## 1. Input Options
With the `deconvolution_spatial` workflow, one or multiple spatial slides can be deconvoluted in one run. For that, a `SpatialData` object for each slide is expected. The spatial slides are deconvoluted **using the same reference**. For the reference, one `MuData` with the gene expression data saved in `mdata.mod["rna"]` is expected as input. Please note, that the same parameter setting is used for each slide. <br> For the **spatial** input, the workflow, therefore, reads in **all `.zarr` objects of a directory** (see below).
<br>

<span class="parameter">input</span><br>
  - <span class="parameter">spatial</span> `String`, Mandatory parameter<br>
        Path to folder containing one or multiple `SpatialDatas` of spatial data. The pipeline is reading in all `SpatialData` files in that folder.

  - <span class="parameter">singlecell</span> `String`, Mandatory parameter<br>
       Path to the MuData **file** (not folder) of the reference single-cell data.


## 2. Cell2Location Options

For each deconvolution method you can specify whether to run it or not: 
<br>

<span class="parameter">run</span> `Boolean`, Default: None<br>
    Whether to run Cell2location


### 2.1 Feature Selection 

You can select genes that are used for deconvolution in two ways. The first option is to provide a reduced feature set as a csv-file that is then used for deconvolution. The second option is to perform gene selection [according to Cell2Location](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html). <br> Please note, that gene selection is **not optional**. If no csv-file is provided, feature selection [according to Cell2Location.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html) is performed. 
<br>

<span class="parameter">feature_selection</span><br>
  - <span class="parameter">gene_list</span> `String`, Default: None<br>
    Path to a csv file containing a reduced feature set. A header in the csv is expected in the first row. All genes of that gene list need to be present in both, spatial slides and scRNA-Seq reference.

  - <span class="parameter">remove_mt</span> `Boolean`, Default: True<br>
    Whether to remove mitochondrial genes from the dataset. This step is performed **before** running gene selection. 

  - <span class="parameter">cell_count_cutoff</span> `Integer`, Default: 15<br>
    All genes detected in less than cell_count_cutoff cells will be excluded. Parameter of the [Cell2Location's gene selection function.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html)

  - <span class="parameter">cell_percentage_cutoff2</span> `Float`, Default: 0.05<br>
    All genes detected in at least this percentage of cells will be included. Parameter of the [Cell2Location's gene selection function.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html)

  - <span class="parameter">nonz_mean_cutoff</span> `Float`, Default: 1.12<br>
    Genes detected in the number of cells between the above-mentioned cutoffs are selected only when their average expression in non-zero cells is above this cutoff.  Parameter of the [Cell2Location's gene selection function.](https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html) 


### 2.2 Reference Model

<span class="parameter">reference</span><br>
  - <span class="parameter">labels_key</span> `String`, Default: None<br>
    Key in `.obs` for label (cell type) information.

  - <span class="parameter">batch_key</span> `String`, Default: None<br>
    Key in `.obs` for batch information.

  - <span class="parameter">layer</span> `String`, Default: None<br>
    Layer in `.layers` to use for the reference model. If None, `.X` will be used. Please note, that Cell2Location expects raw counts as input.

  - <span class="parameter">categorical_covariate_key</span> `String`, Default: None<br>
    Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space).

  - <span class="parameter">continuous_covariate_keys</span> `String`, Default: None<br>
    Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space)

  - <span class="parameter">max_epochs</span> `Integer`, Default: _np.min([round((20000 / n_cells) * 400), 400])_<br>
    Number of epochs.

  - <span class="parameter">use_gpu</span> `Boolean`, Default: True<br>
    Whether to use GPU for training.
   


### 2.3 Spatial Model 

<span class="parameter">spatial</span><br>
  - <span class="parameter">batch_key</span> `String`, Default: None<br>
    Key in `.obs` for batch information.

  - <span class="parameter">layer</span> `String`, Default: None<br>
    Layer in `.layers` to use for the reference model. If None, `.X` will be used. Please note, that Cell2Location expects raw counts as input.

  - <span class="parameter">categorical_covariate_key</span> `String`, Default: None<br>
    Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space).

  - <span class="parameter">continuous_covariate_keys</span> `String`, Default: None<br>
    Comma-separated without spaces, e.g. _key1,key2,key3_. Keys in `.obs` that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space)

  - <span class="parameter">N_cells_per_location</span> `Integer`, Mandatory parameter<br>
    Expected cell abundance per voxel. Please refer to the [Cell2Location documentation](https://cell2location.readthedocs.io/en/latest/index.html) for more information. 

  - <span class="parameter">detection_alpha</span> `Float`, Mandatory parameter<br>
    Regularization of with-in experiment variation in RNA detection sensitivity. Please refer to the [Cell2Location documentation](https://cell2location.readthedocs.io/en/latest/index.html) for more information. 

  - <span class="parameter">max_epochs</span> `Integer`, Default: _np.min([round((20000 / n_cells) * 400), 400])_<br>
    Number of epochs.

  - <span class="parameter">use_gpu</span> `Boolean`, Default: True<br>
    Whether to use GPU for training.


<br>


You can specify whether both models (spatial and reference) should be saved with the following parameter:
<br>

<span class="parameter">save_models</span>, Default: False<br>
    Whether to save the reference & spatial mapping models.

<span class="parameter">export_gene_by_spot</span>, Default: False<br>
    Whether to save a gene by spot matrix for each cell type in a layer.


## 3. Tangram Options

For each deconvolution method you can specify whether to run it or not: 
<br>

<span class="parameter">run</span> `Boolean`, Default: None<br>
    Whether to run Tangram


### 3.1 Feature Selection

You can select genes that are used for deconvolution in two ways. The first option is to provide a reduced feature set as a csv-file that is then used for deconvolution. The second option is to perform gene selection via [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) **on the reference scRNA-Seq data**, as [suggested by Tangram](https://tangram-sc.readthedocs.io/en/latest/tutorial_sq_link.html#Pre-processing). The top `n_genes` of each group make up the reduced gene set. <br> Please note, that gene selection is **not optional**. If no csv-file is provided, feature selection via [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html) is performed. 
<br>

<span class="parameter">feature_selection</span><br>
  - <span class="parameter">gene_list</span> `String`, Default: None<br>
    Path to a csv file containing a reduced feature set. A header in the csv is expected in the first row. All genes of that gene list need to be present in both, spatial slides and scRNA-Seq reference.

___Parameters for `scanpy.tl.rank_genes_groups` gene selection___ 
  - <span class="parameter">rank_genes</span> <br>
    - <span class="parameter">labels_key</span> `String`, Default: None<br>
        Which column in `.obs` of the reference to use for the `groupby` parameter of [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html).

    - <span class="parameter">layer</span> `String`, Default: None<br>
        Which layer of the reference to use for [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html). If None, `.X` is used.

    - <span class="parameter">n_genes</span> `Integer`, Default: 100<br>
        How many top genes to select of each `groupby` group.

    - <span class="parameter">test_method</span>  `['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']`, Default: 't-test_overestim_var'<br>
        Which test method to use.

    - <span class="parameter">correction_method</span> `['benjamini-hochberg', 'bonferroni']`, Default: ' benjamini-hochberg'<br>
        Which p-value correction method to use. Used only for 't-test', 't-test_overestim_var', and 'wilcoxon'. 


### 3.2 Model

<span class="parameter">model</span><br>
  - <span class="parameter">labels_key</span> `String`, Default: None<br>
    Key in `.obs` for label (cell type) information.

  - <span class="parameter">num_epochs</span> `Integer`, Default: 1000<br>
    Number of epochs.

  - <span class="parameter">device</span> `String`, Default: 'cpu'<br>
    Which device to use.

  - <span class="parameter">kwargs</span><br>
    In `kwargs`, the user has the possibility to specify parameters for [tangram.mapping_utils.map_cells_to_space](https://tangram-sc.readthedocs.io/en/latest/classes/tangram.mapping_utils.map_cells_to_space.html?highlight=mapping_utils%20map_cells_to_space#tangram.mapping_utils.map_cells_to_space). You can add or remove any parameters of the function.
   

