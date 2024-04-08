<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# Visualization YAML

In this documentation, the parameters of the `visualization` configuration yaml file are explained. 
This file is generated running `panpipes vis config`.  <br> The individual steps run by the pipeline are described in the [visualization workflow](WHEERREEEE)??????
When running the visualization workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.
However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)

You can download the different ingestion `pipeline.yml` files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes vis config`: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_vis/pipeline.yml)
- `pipeline.yml` file for [Visualizing data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/visualization/pipeline_yml.html): [Download here](https://panpipes-tutorials.readthedocs.io/en/latest/_downloads/29daa86241829b362152785caf30ab61/pipeline.yml)

## Compute resources options
<span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.
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

<span class="parameter">condaenv</span> `String` (Path)<br>
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment

## Loading and merging data options
### Data format

<span class="parameter">sample_prefix</span> `String`, Mandatory parameter, Default: test<br>
Prefix for the sample that comes out of the filtering/ preprocessing steps of the workflow.

<span class="parameter">mudata_obj</span> `String`, Mandatory parameter <br>
 Path to the output file from preprocessing (e.g. `../vis/test.h5mu`).
 Ensure that the submission file is in the right format and that the correct path is provided.

<span class="parameter">modalities</span><br>
<span class="parameter">rna</span> `Boolean`, Default: True <br>
<span class="parameter">prot</span> `Boolean`, Default: True <br>
<span class="parameter">atac</span> `Boolean`, Default: False <br>
<span class="parameter">rep</span> `Boolean`, Default: True <br>
<span class="parameter">multimodal</span> `Boolean`, Default: True <br>
Set the modalities to True or False depending on what is present in the mudata_obj

<span class="parameter">grouping_vars</span> `String`, Default: sample_id  rna:leiden_res0.6 <br>
On dot plots and bar plots, grouping vars are used to group other features (for categorical, continuous, and feature plots).

## Plot Markers 
<span class="parameter">custom_markers</span><br>
  - <span class="parameter">files</span><br>
  
The csv files for full and minimal must contain three columns:

  | mod  | feature  | group        |
  |------|----------|--------------|
  | prot | prot_CD8 | Tcellmarkers |
  | rna  | CD8A     | Tcellmarkers |

  - <span class="parameter">full:</span><br>
The full list will be plotted in dot plots and matrix plots, with one plot per group. 

 - <span class="parameter">minimal:</span><br>
The shorter list will be plotted on umaps as well as other plot types, with one plot per group. 

 | feature_1 | feature_2 | colour         |
 |-----------|-----------|----------------|
 | CD8A      | prot_CD8  |                |
 | CD4       | CD8A      | doublet_scores |
    
  
- <span class="parameter">paired_scatter:</span>`String`, Default:  <br>
  Where different normalization exists in a modality, choose which one to use, set X or leave blank to use the mdata[mod].X assay. 

- <span class="parameter">layers:</span><br>
  - <span class="parameter">rna:</span>`String`, Default: logged_counts<br>
  - <span class="parameter">prot:</span>`String`, Default: logged_counts<br> CHEEECKKKK
  - <span class="parameter">atac:</span>`String`, Default: logged_counts<br> CHHHHECKKK
    
## Plot metadata variables 

- <span class="parameter">categorical_vars:</span>`String`, Default: &categorical_vars<br>
  - <span class="parameter">all:</span>`String`, Default: rep:receptor_subtype  sample_id<br>
Metrics to be plotted on every modality. 
  - <span class="parameter">rna:</span>`String`, Default: rna:predicted_doublets  rna:phase<br>
  - <span class="parameter">prot:</span>`String`, Default: prot:leiden_res0.2    prot:leiden_res1<br>
  - <span class="parameter">atac:</span>`String`, Default: <br>
  - <span class="parameter">rep:</span>`String`, Default: prot:rep:has_ir<br>
  - <span class="parameter">multimodal:</span>`String`, Default: leiden_totalVI    mdata_colsr<br>

- <span class="parameter">continuous_vars:</span>`String`, Default: &continuous_vars<br>
  - <span class="parameter">all:</span>`String`, Default:leiden_res0.5<br>
Metrics to be plotted on every modality. 
  - <span class="parameter">rna:</span>`String`, Default: rna:total_counts<br>
  - <span class="parameter">prot:</span>`String`, Default: prot:total_counts<br>
  - <span class="parameter">atac:</span>`String`, Default: <br>
  - <span class="parameter">multimodal:</span>`String`, Default: rna:total_counts    prot:total_counts<br>
  
- <span class="parameter"paired_scatter:</span>`String`, Default: scatter_features.csv<br>
The scatter_features.csv file should have the following format:

 | feature_1 | feature_2 | colour         |
 |-----------|-----------|----------------|
 |rna:total_counts | prot:total_counts  | doublet_scores|

## Plot style 


 
      



