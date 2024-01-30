
# Spatial QC YAML 

In this documentation, the parameters of the `qc_spatial` yaml file are explained. In general, the user can leave parameters empty to use defaults. <br>  The individual steps run by the pipeline are described in the [spatial QC workflow](../workflows/ingest_spatial.md). 



## 0. Compute Resource Options

| `resources` |  |
| --- | --- |
| `threads_high` | __`int`__ (default: 1) <br> Number of threads used for high intensity computing tasks. For each thread, there must be enough memory to load all input files and create MuDatas.  |
| `threads_medium` | __`int`__ (default: 1) <br> Number of threads used for medium intensity computing tasks. For each thread, there must be enough memory to load your mudata and do computationally light tasks. |
| `threads_low` | __`int`__ (default: 1) <br> Number of threads used for low intensity computing tasks. For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.|



|  |  |
| ---- | --- |
| `condaenv` | __`str`__ (default: None) <br> Path to conda environment that should be used to run `Panpipes`. Leave blank if running native or your cluster automatically inherits the login node environment |




## 1. Loading Options 

|  |  |
| ---- | --- |
| `project` | __`str`__ (default: None) <br> Project name |
| `submission_file` | __`str`__ (not optional) <br> Path to the submission file. The submission file specifies the input files. Please refer to the [general guidelines](../usage/setup_for_spatial_workflows.md) for details on the format of the file. |


## 2. QC Options 
This part of the workflow allows to generate additional QC metrics that can be used for filtering/preprocessing. Basic QC metrics using [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html) are in every case calculated. The computation of the additional QC metrics is **optional**. Please, leave the parameters empty to avoid running.

|  |  |
| ---- | --- |
| `ccgenes` | __`str`__ (default: None) <br> Path to tsv-file used to run the function [scanpy.tl.score_genes_cell_cycle](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html). It is expected, that the tsv-file has two columns with names `cc_phase` and `gene_name`. `cc_phase` can either be `s` or `g2m`.**Varying the column names or `cc_phase` values will result in an error.** Please refer to the [general guidelines](../usage/gene_list_format.md) for more information on the tsv file. <br> Instead of a path, the user can specify the parameter as "default" which then uses a [provided tsv file](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/cell_cycle_genes.tsv).|
| `custom_genes_file` | __`str`__ (default: None) <br> Path to csv-file containing a gene list with columns `group` and `feature`. **Varying the column names will result in an error.** Please refer to the [general guidelines](../usage/gene_list_format.md) for more information about the file. <br> The gene list is used to calculate the proportions of genes of a group in the cells/spots. More precise, the groups & genes are used for the `qc_vars` parameter of the function [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html) which accordingly calculates proportions. <br> Additionally the gene list is used to compute gene scores with [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) <br> Instead of a path, the user can specify the parameter as "default" which then uses a [provided csv file](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/qc_genelist_1.0.csv).|
| `calc_proportions` | __`str`__ (default: None) <br>  Comma-separated string without spaces, e.g. _mito,hp,rp_. <br> For which groups of the csv-file specified in `custom_genes_file` to calculate percentages. |
| `score_genes` | __`str`__ (default: None) <br> Comma-separated string without spaces, e.g. _mito,hp,rp_. <br> For which groups of the csv-file specified in `custom_genes_file`  to run [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) |

The following parameters specify the QC metrics to plot in violin and spatial embedding plots. Plots are generated for each slide specified in the submission file separately. 

| `plotqc` |  |
| --- | --- |
| `grouping_var` | __`str`__ (default: None) <br> Comma-separated string without spaces, e.g. _sample_id,batch_ of categorical columns in `.obs`. One violin will be created for each group in the violin plot. Not mandatory, can be left empty.|
| `spatial_metrics` | __`str`__ (default: None) <br>  Comma-separated string without spaces, e.g. _total_counts,n_genes_by_counts_ of columns in `.obs` or `.var`. <br>Specifies which metrics to plot. If metric is present in both, `.obs` and `.var`, **both will be plotted.**|
