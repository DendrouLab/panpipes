<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style> 

# Refmap YAML 
In this documentation, the parameters of the `refmap` configuration yaml file are explained. 
This file is generated running `panpipes refmap config`.  <br> The individual steps run by the pipeline are described in the [Reference Mapping workflow](https://github.com/DendrouLab/panpipes/blob/main/docs/workflows/refmap.md). 


When running the refmap workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.
However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-tutorials.readthedocs.io/en/latest/refmap_pancreas/Reference_mapping.html)

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)

You can download the different refmap `pipeline.yml` files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes refmap config: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_refmap/pipeline.yml)
- `pipeline.yml` file for [Reference Mapping Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/refmap_pancreas/Reference_mapping.html): [Download here](https://panpipes-tutorials.readthedocs.io/en/latest/_downloads/cfb2a3d64a5e7b2cabe7ee8e1ac5fe61/pipeline.yml)


## Compute resources options

<span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.  Check [threads_tasks_panpipes.csv](https://github.com/DendrouLab/panpipes/blob/threads_doc_g/docs/yaml_docs/threads_tasks_panpipes.csv) for more information on which threads each specific task requires.
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

  - <span class="parameter">condaenv</span> `String` (Path)<br>
Path to conda environment that should be used to run panpipes.

  - <span class="parameter">queues:</span> `String` (Path)<br>
In case a special queue is required for long jobs or if the user has access to a GPU-specific queue. Otherwise, leave it blank. 
    - <span class="parameter">long:</span> `String` (Path)<br>
    - <span class="parameter">gpu:</span> `String` (Path)<br>


## Loading data options
### Query Dataset

- <span class="parameter">query</span> `String`, Default: path/to/data<br>
    Give the path to the desired data. Formats accepted include raw10x, preprocessed quality filtered mudata or anndata as input query
- <span class="parameter">modality</span> `String`, Default: rna<br>
If mudata was provided then specify the modality to be used. Currently, only RNA modality is supported. 
- <span class="parameter">query_batch</span> `String`, Default: <br>
Only to be filled if the data provided had a batch correction, if so specify the column this is in. If not, leave blank 
- <span class="parameter">query_celltype</span> `String`, Default: <br>
If the query provided has celltype annotations that should be compared to the transferred labels. If not, leave blank.

## Scvi tools parameters 

- <span class="parameter">reference_data</span> `String`, Default: path/to/mudata<br>
Specify one or more reference models to be used as reference. Users can also specify their own reference built using `pipeline_integration`.
Leave blank for no model specification.

- <span class="parameter">totalvi:</span> `String`, Default: path/to/totalvi<br>
Provide path to totalvi saved model. Multiple paths can be provided as a list:
```yaml 
totalvi: 
  - path_to_totalvi1
  - path_to_totalvi2

```
  -

- <span class="parameter">impute_proteins</span> `Boolean`, Default: False<br> 
- <span class="parameter">transform_batch</span> `String`, Default:<br>
Transform_batch is a batch-covariate specific to totalvi, allows the model to use the batch information in the query to mitigate 
differences in protein sequencing depth.
- <span class="parameter">scvi</span> `String`, Default: path/to/scvi Mandatory, Provide a path to the scvi model. Multiple paths can be provided as a list: <br>

```yaml 
scvi: 
  - path_to_totalvi1
  - path_to_totalvi2

```

- <span class="parameter">scanvi</span> `String`, Default:path/to/scanvi Mandatory, Provide a path to the scvi model.<br>
- <span class="parameter">run_randomforest</span> `Boolean`, Default:False<br>
Set to true if the reference model has a trained random forest classifier to transfer the labels. 

## Training parameters 
To reuse the same params in multiple locations, please use anchors (&) and scalars (*) in the relevant place, i.e. if specifying &rna_neighbors, the same params will be called by *rna_neighbors where referenced. Check our documentation for more info on using anchors and scalars

- <span class="parameter">training_plan:</span><br>
  - <span class="parameter">totalvi:</span> Default: array of training parameters. <br>For the full list of parameters check [here](https://docs.scvi-tools.org/en/0.14.1/api/reference/scvi.model.TOTALVI.train.html). to reuse the same parameters in other locations use an anchor, for example writing `totalvi: &totalvitraining` and will ensure the same array is reused when referencing it as `*totalvitraining`. In this example the `&totalvitraining` array contains the two parameters `max_epochs` and `weight_decay` 
    - <span class="parameter">max_epochs</span> `Integer`, Default: 200<br>
    - <span class="parameter">weight_decay</span> `Float`, Default: 0.0<br>
    Recommended weight decay is 0.0. This ensures the latent representation of the reference cells will remain exactly the same if passing them through this new query model.
    - <span class="parameter">scvi</span> Array of training parameters, Default: `*totalvitraining` (reuse the same array as specified above)<br>
    - <span class="parameter">scanvi</span> Array of training parameters, Default: `*totalvitraining` (reuse the same array as specified above) <br>

## Neighbors parameters to calculate umaps 
This can be on either query alone, or query+ reference dataset. 

- <span class="parameter">neighbors:</span><br>
  - <span class="parameter">npcs</span> `Integer`, Default: 30<br>
Number of Principal Components to calculate for neighbours and umap. If no correction is applied, PCA will be calculated and used to run UMAP and clustering on.
And if Harmony is the method of choice, it will use these components to create a corrected dim red.
  - <span class="parameter">k</span> `Integer`, Default: 30<br>
This is the number of neighbours
  - <span class="parameter">metric</span> `String`, Default: euclidean<br>
Options here include cosine and euclidean
  - <span class="parameter">method</span> `String`, Default: sanpy<br>
Options here include scanpy, and hnsw (from scvelo)

## Run scib metrics on query
Running scib on query data after transferring labels, where available (with the totalvi and scanvi models), or using default leiden clustering after training the vae model (scvi)
Check [documentation](https://scib.readthedocs.io/en/latest/) for the metrics used 
- <span class="parameter">scib:</span><br>
  - <span class="parameter">run</span> `Boolean`, Default: False<br>
  - <span class="parameter">cluster_key</span> `String`, Default: predictions<br>
Used for ARI and NMI, if left empty will default to leiden clustering calculated on the new latent representation after reference mapping.
  - <span class="parameter">batch_key</span> `String`, Default: <br>
 Used for clisi_graph_embed and if no batch is present the metrics will not be included in the results. If left blank will default do cluster_key defauls.
  - <span class="parameter">celltype_key</span> `String`, Default: celltype <br>






