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

