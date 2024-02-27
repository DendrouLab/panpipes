<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# Clustering YAML 

In this documentation, the parameters of the `clustering` configuration yaml file are explained.
This file is generated running `panpipes clustering config`. <br>
The individual steps run by the pipeline are described in [clustering worlfow](https://panpipes-pipelines.readthedocs.io/en/latest/workflows/clustering.html)

When running the clustering workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.

However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).
For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)

You can download the different clustering pipeline.yml files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes clustering config`: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_clustering/pipeline.yml)
- `pipeline.yml` for [Clustering Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/_downloads/3895aa0ba60017b15ee1aa6531dc8c25/pipeline.ym)

## Compute resources options

- <span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.
Specified by the following three parameters:
  - <span class="parameter">threads_high</span> `Integer`, Default: 2<br>
        Number of threads used for high intensity computing tasks. 
        For each thread, there must be enough memory to load all your input files at once and create the MuData object.

  - <span class="parameter">threads_medium</span> `Integer`, Default: 2<br>
        Number of threads used for medium intensity computing tasks.
        For each thread, there must be enough memory to load your mudata and do computationally light tasks.

  - <span class="parameter">threads_low</span> `Integer`, Default: 2<br>
  	    Number of threads used for low intensity computing tasks.
        For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.
  - <span class="parameter">fewer_jobs</span> `Boolean`, Default: True<br>
  
  - <span class="parameter">condaenv</span> `String` (Path)<br>
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment
