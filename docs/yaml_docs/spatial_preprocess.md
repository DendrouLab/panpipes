<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# Spatial Preprocessing YAML

In this documentation, the parameters of the `preprocess_spatial` configuration yaml file are explained. 
This file is generated running `panpipes preprocess_spatial config`.  <br> The individual steps run by the pipeline are described in the [spatial preprocessing workflow](../workflows/preprocess_spatial.md). 


For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)

When running the preprocess workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.
However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).
You can download the different preprocess pipeline.yml files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes preprocess_spatial config`: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_preprocess_spatial/pipeline.yml)
- `pipeline.yml` file for [Preprocessing spatial data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/preprocess_spatial_data/preprocess_spatial_data_with_panpipes.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/preprocess_spatial_data/pipeline.yml)



## 0. Compute Resource Options

<span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.  Check [threads_tasks_panpipes](./threads_tasks_panpipes.md) for more information on which threads each specific task requires.
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

With the preprocess_spatial workflow, one or multiple `SpatialData` objects can be preprocessed in one run. The workflow **reads in all `.zarr` objects of a directory**. The `SpatialData` objects in the directory need to be of the same assay (Vizgen, Visium, or Xenium). The workflow then runs the preprocessing of each `SpatialData` object separately with the same parameters that are specified in the yaml file. 
<br>

<span class="parameter">input_dir</span> `String`, Mandatory parameter<br>
    Path to the folder containing all input `zarr` files.

<span class="parameter">assay</span> [`'visium'`, `'vizgen'`], Default: `'visium'`<br>
     Spatial transcriptomics assay of the `zarr` files in `input_dir`.



## 2. Filtering Options

<span class="parameter">filtering</span><br>
  - <span class="parameter">run</span> `Boolean`, Default: False<br>
        Whether to run filtering. **If `False`, will not filter the data and will not produce post-filtering plots.**

  - <span class="parameter">keep_barcodes</span> `String`, Default: None<br>
        Path to a csv-file that has **no header** containing barcodes you want to keep. Barcodes that are not in the file, will be removed from the dataset before filtering the dataset with the thresholds specified below. 
<br>


With the parameters below you can specify thresholds for filtering. The filtering is fully customisable to any columns in `.obs` or `.var`. You are not restricted by the columns given as default. When specifying a column name, please make sure it exactly matches the column name in the table of the `SpatialData` object. <br> Please also make sure, that the specified metrics are present in all `SpatialData` objects of the `input_dir`, i.e. the `SpatialData` objects for that the preprocessing is run.


---
    spatial:
        obs:
            min:
                total_counts: 
            max:
                pct_counts_mt:
            bool: 
        var:
            min:
                n_cells_by_counts: 
            max:
                total_counts:
---


## 3. Post-Filter Plotting

The parameters below specify which metrics of the filtered data to plot. As for the [QC](./spatial_qc.md), violin and spatial embedding plots are generated for each slide separately. 
<br>

<span class="parameter">plotqc</span><br>
  - <span class="parameter">grouping_var</span> `String`, Default: None<br>
        Comma-separated string without spaces, e.g. _sample_id,batch_ of categorical columns in `.obs`. One violin will be created for each group in the violin plot. Not mandatory, can be left empty.

  - <span class="parameter">spatial_metrics</span> `String`, Default: None<br>
        Comma-separated string without spaces, e.g. _total_counts,n_genes_by_counts_ of columns in `.obs` or `.var`. <br>Specifies which metrics to plot. If metric is present in both, `.obs` and `.var`, **both will be plotted.**
    

## 4. Normalization, HVG Selection, and PCA Options

### **4.1 Normalization and HVG Selection**
`Panpipes` offers two different normalization and HVG selection flavours, `'seurat'` and `'squidpy'`. <br> The `'seurat'`  flavour first selects HVGs on the raw counts using analytic Pearson residuals, i.e. [scanpy.experimental.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.highly_variable_genes.html). Afterwards, analytic Pearson residual normalization is applied, i.e. [scanpy.experimental.pp.normalize_pearson_residuals](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html). Parameters of both functions can be specified by the user in the yaml file. <br>The `'squidpy'` flavour runs the basic scanpy normalization and HVG selection functions, i.e. [scanpy.pp.normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html), [scanpy.pp.log1p](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html), and [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html). 
<br> 

<span class="parameter">norm_hvg_flavour</span>[`'squidpy'`, `'seurat'`], Default: None<br>
    Normalization and HVG selection flavour to use. If None, will not run normalization nor HVG selection. 
<br>

___Parameters for `norm_hvg_flavour` == `'squidpy'`___ <br>

<span class="parameter">squidpy_hvg_flavour</span>[`'seurat'`,`'cellranger'`,`'seurat_v3'`], Default: 'seurat'<br>
    Flavour to select HVGs, i.e.`flavor` parameter of the function [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html).

<span class="parameter">min_mean</span>`Float`, Default: 0.05<br>
    Parameter in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html).

<span class="parameter">max_mean</span>`Float`, Default: 1.5<br>
    Parameter in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html). 

<span class="parameter">min_disp</span>`Float`, Default: 0.5<br>
    Parameter in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html).

___Parameters for `norm_hvg_flavour` == `'seurat'`___ <br>

<span class="parameter">theta</span>`Float`, Default: 100<br>
    The negative binomial overdispersion parameter for pearson residuals. The same value is used for [HVG selection]((https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.highly_variable_genes.html)) and [normalization](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html). 

<span class="parameter">clip</span>`Float`, Default: None<br>
    Specifies clipping of the residuals. <br>`clip` can be specified as: <br> <ul><li> <u>None</u>: residuals are clipped to the interval [-sqrt(n_obs), sqrt(n_obs)] </li><li><u>A float value</u>: if float c specified: clipped to the interval [-c, c]</li> <li> <u>np.Inf</u>: no clipping</li></ul> 

___Parameters for both `norm_hvg_flavour` flavours___ <br>

<span class="parameter">n_top_genes</span>`Integer`, Default: 2000<br>
    Number of genes to select. Mandatory for `norm_hvg_flavour='seurat'` and `squidpy_hvg_flavour='seurat_v3'`.

<span class="parameter">filter_by_hvg</span>`Boolean`, Default: False<br>
    Subset the data to the HVGs. 

<span class="parameter">hvg_batch_key</span>`String`, Default: None<br>
    If specified, HVGs are selected within each batch separately and merged. 


### **4.2 PCA**

After normalization and HVG selection, PCA is run and the PCA and elbow plot are plotted. For that, the user can specify the number of PCs for the PCA computation and for the elbow plot, i.e. the same number is used for both. 
<br>

<span class="parameter">n_pcs</span>`Integer`, Default: 50<br>
    Number of PCs to compute.

## 5. Concatenation

In case multiple `SpatialData` objects have been preprocessed separately, the user has the option to concatenate the preprocessed objects in the end. 
<br>

<span class="parameter">concat</span>`Boolean`, Default: False<br>
    Whether to concatenate all preprocessed `SpatialData` objects. The concatenated object is saved to `./concatenated.data/concatenated.zarr`




