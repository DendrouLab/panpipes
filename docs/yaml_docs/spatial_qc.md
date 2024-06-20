<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# Spatial QC YAML 

In this documentation, the parameters of the `qc_spatial` configuration yaml file are explained. 
This file is generated running `panpipes qc_spatial config`.  <br> The individual steps run by the pipeline are described in the [spatial QC workflow](../workflows/ingest_spatial.md). 

When running the qc workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.
However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).
You can download the different ingestion pipeline.yml files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes qc_spatial config`: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_qc_spatial/pipeline.yml)
- `pipeline.yml` file for [Ingesting 10X Visium data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_visium_data/Ingesting_visium_data_with_panpipes.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/ingesting_visium_data/pipeline.yml)
- `pipeline.yml` file for [Ingesting MERFISH data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_merfish_data/Ingesting_merfish_data_with_panpipes.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/ingesting_merfish_data/pipeline.yml)

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)


## 0. Compute Resource Options


<span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.  check [threads_tasks_panpipes.csv](https://github.com/DendrouLab/panpipes/blob/threads_doc_g/docs/yaml_docs/threads_tasks_panpipes.csv) for more information on which threads each specific task requires.
Specified by the following three parameters:
  - <span class="parameter">threads_high</span> `Integer`, Default: 1<br>
        Number of threads used for high intensity computing tasks. 
        For each thread, there must be enough memory to load all your input files at once and create the MuData object.

  - <span class="parameter">threads_medium</span> `Integer`, Default: 1<br>
        Number of threads used for medium intensity computing tasks.
        For each thread, there must be enough memory to load your mudata and do computationally light tasks.

  - <span class="parameter">threads_low</span> `Integer`, Default: 1<br>
  	    Number of threads used for low intensity computing tasks.
        For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.

<span class="parameter">condaenv</span> `String`<br>
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment




## 1. Loading Options 

<span class="parameter">project</span> `String`, Default: None<br>
    Project name.

<span class="parameter">submission_file</span> `String`, Mandatory parameter<br>
   Path to the submission file. The submission file specifies the input files. Please refer to the [general guidelines](../usage/setup_for_spatial_workflows.md) for details on the format of the file.


## 2. QC Options 
This part of the workflow allows to generate additional QC metrics that can be used for filtering/preprocessing. Basic QC metrics using [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html) are in every case calculated. The computation of the additional QC metrics is **optional**. Please, leave the parameters empty to avoid running.
<br>

<span class="parameter">ccgenes</span> `String`, Default: None<br>
    Path to tsv-file used to run the function [scanpy.tl.score_genes_cell_cycle](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html). It is expected, that the tsv-file has two columns with names `cc_phase` and `gene_name`. `cc_phase` can either be `s` or `g2m`.**Varying the column names or `cc_phase` values will result in an error.** Please refer to the [general guidelines](../usage/gene_list_format.md) for more information on the tsv file. <br> Instead of a path, the user can specify the parameter as "default" which then uses a [provided tsv file](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/cell_cycle_genes.tsv).

<span class="parameter">custom_genes_file</span> `String`, Default: None<br>
     Path to csv-file containing a gene list with columns `group` and `feature`. **Varying the column names will result in an error.** Please refer to the [general guidelines](../usage/gene_list_format.md) for more information about the file. <br> The gene list is used to calculate the proportions of genes of a group in the cells/spots. More precise, the groups & genes are used for the `qc_vars` parameter of the function [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html) which accordingly calculates proportions. <br> Additionally the gene list is used to compute gene scores with [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) <br> Instead of a path, the user can specify the parameter as "default" which then uses a [provided csv file](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/qc_genelist_1.0.csv).

<span class="parameter">calc_proportions</span> `String`, Default: None<br>
     Comma-separated string without spaces, e.g. _mito,hp,rp_. <br> For which groups of the csv-file specified in `custom_genes_file` to calculate percentages. 

<span class="parameter">score_genes</span> `String`, Default: None<br>
    Comma-separated string without spaces, e.g. _mito,hp,rp_. <br> For which groups of the csv-file specified in `custom_genes_file`  to run [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) 

<br>
The following parameters specify the QC metrics to plot in violin and spatial embedding plots. Plots are generated for each slide specified in the submission file separately. <br>
<br>

<span class="parameter">plotqc</span><br>
  - <span class="parameter">grouping_var</span> `String`, Default: None<br>
        Comma-separated string without spaces, e.g. _sample_id,batch_ of categorical columns in `.obs`. One violin will be created for each group in the violin plot. Not mandatory, can be left empty.

  - <span class="parameter">spatial_metrics</span> `String`, Default: None<br>
        Comma-separated string without spaces, e.g. _total_counts,n_genes_by_counts_ of columns in `.obs` or `.var`. <br>Specifies which metrics to plot. If metric is present in both, `.obs` and `.var`, **both will be plotted.**

