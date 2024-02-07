
# Spatial Preprocessing YAML

In this documentation, the parameters of the `preprocess_spatial` yaml file are explained. In general, the user can leave parameters empty to use defaults. <br> The individual steps run by the pipeline are described in the [spatial preprocess workflow](../workflows/preprocess_spatial.md). 



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

With the preprocess_spatial workflow, one or multiple `MuData` objects can be preprocessed in one run. The workflow **reads in all `.h5mu` objects of a directory**. The `MuData` objects in the directory need to be of the same assay (vizgen or visium). The workflow then runs the preprocessing of each `MuData` object separately with the same parameters that are specified in the yaml file. 

| |  |
| ---- | --- |
| `input_dir` | __`str`__ (not optional) <br> Path to the folder containing all input `h5mu` files. |
| `assay` | __`str` [`'visium'`, `'vizgen'`]__ (default: 'visium') <br> Spatial transcriptomics assay of the `h5mu` files in `input_dir`.|


## 2. Filtering Options


| `filtering` |  |
| --- | --- | 
| `run` | __`bool`__ (default: False) <br> Whether to run filtering. **If `False`, will not filter the data and will not produce post-filtering plots.** |  
| `keep_barcodes` | __`str`__ (default: None) <br> Path to a csv-file that has **no header** containing barcodes you want to keep. Barcodes that are not in the file, will be removed from the dataset before filtering the dataset with the thresholds specified below. |


With the parameters below you can specify thresholds for filtering. The filtering is fully customisable to any columns in `.obs` or `.var`. You are not restricted by the columns given as default. When specifying a column name, please make sure it exactly matches the column name in the h5mu object. <br> Please slso make sure, that the specified metrics are present in all `h5mu` objects of the `input_dir`, i.e. the `MuData` objects for that the preprocessing is run.


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

| `plotqc` |  |
| --- | --- | 
| `grouping_var` | __`str`__ (default: None) <br>  Comma-separated string without spaces, e.g. _sample_id,batch_ of categorical columns in `.obs`. One violin will be created for each group in the violin plot. Not mandatory, can be left empty. |  
| `spatial_metrics` | __`str`__ (default: None) <br>  Comma-separated string without spaces, e.g. _total_counts,n_genes_by_counts_ of columns in `.obs` or `.var`. <br>Specifies which metrics to plot. If metric is present in both, `.obs` and `.var`, **both will be plotted.** |


## 4. Normalization, HVG Selection, and PCA Options

### **Normalization and HVG Selection** <br>

`Panpipes` offers two different normalization and HVG selection flavours, `'seurat'` and `'squidpy'`. <br> The `'seurat'`  flavour first selects HVGs on the raw counts using analytic Pearson residuals, i.e. [scanpy.experimental.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.highly_variable_genes.html). Afterwards, analytic Pearson residual normalization is applied, i.e. [scanpy.experimental.pp.normalize_pearson_residuals](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html). Parameters of both functions can be specified by the user in the yaml file. <br>The `'squidpy'` flavour runs the basic scanpy normalization and HVG selection functions, i.e. [scanpy.pp.normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html), [scanpy.pp.log1p](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html), and [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html). <br> 


| `spatial` |  |
| --- | --- |
| `norm_hvg_flavour` | __`str` [`'squidpy'`, `'seurat'`]__ (default: None) <br>  Normalization and HVG selection flavour to use. If None, will not run normalization nor HVG selection. | 

___Parameters for `norm_hvg_flavour` == `'squidpy'`___ 
| |  |
| --- | --- |
| `squidpy_hvg_flavour` | __`str` [`'seurat'`,`'cellranger'`,`'seurat_v3'`]__ (default: 'seurat') <br>   Flavour to select HVGs, i.e.`flavor` parameter of the function [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html).|
| `min_mean` | __`float`__ (default: 0.05) <br>  Parameter in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html).|
| `max_mean` | __`float`__ (default: 1.5) <br>  Parameter in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html). |
| `min_disp` | __`float`__ (default: 0.5) <br>  Parameter in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html). |

___Parameters for `norm_hvg_flavour` == `'seurat'`___ 
| |  |
| --- | --- |
| `theta` | __`float`__ (default: 100) <br>  The negative binomial overdispersion parameter for pearson residuals. The same value is used for [HVG selection]((https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.highly_variable_genes.html)) and [normalization](https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html).  |
| `clip` | __`float`__ (default: None) <br> Specifies clipping of the residuals. <br>`clip` can be specified as: <br> <ul><li> <u>None</u>: residuals are clipped to the interval [-sqrt(n_obs), sqrt(n_obs)] </li><li><u>A float value</u>: if float c specified: clipped to the interval [-c, c]</li> <li> <u>np.Inf</u>: no clipping</li></ul> | 

___Parameters for both `norm_hvg_flavour` flavours___ 
| |  |
| --- | --- |
| `n_top_genes` | __`int`__ (default: 2000) <br> Number of genes to select. Mandatory for `norm_hvg_flavour='seurat'` and `squidpy_hvg_flavour='seurat_v3'`.|
| `filter_by_hvg` | __`bool`__ (default: False) <br> Subset the data to the HVGs. <br>  |
| `hvg_batch_key` | __`str`__ (default: None) <br> If specified, HVGs are selected within each batch separately and merged. |


### **PCA**

After normalization and HVG selection, PCA is run and the PCA and elbow plot are plotted. For that, the user can specify the number of PCs for the PCA computation and for the elbow plot, i.e. the same number is used for both. 

| |  |
| --- | --- |
| `n_pcs` | __`int`__ (default: 50) <br> Number of PCs to compute. |

