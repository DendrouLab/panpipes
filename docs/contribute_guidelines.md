# Contributing to panpipes

### 1. Explaining the Magic

`panpipes` is written in python using [`cgat-core`](https://github.com/cgat-developers/cgat-core) a workflow management system that allows users to create and execute complex data analysis pipelines with ease.
The pipeline is designed in modular fashion, with its [9 workflows](../README.md) for the comprehensive analysis of single cell multiomics and spatial transcriptomics data.

Each step of the analysis is performed by an independent python script importing the relevant libraries. Wherever possible, `panpipes` parallelizes the computation of multiple tasks, such as clustering and marker discovery on a single cell experiment at multiple resolutions, or performing batch correction and multimodal integration methods simultaneously (where applicable). 

Each of the workflows requires a configuration yml file `pipeline.yml` which can be generated running `panpipes NAME_OF_WORKFLOW config`. This file can be customized to taylor the needs and questions of the analyst and it's required to initiate and configure the workflows appropriately. (See the tutorials for more information)

Let's check for example the batch_correction method harmony, which we include in the script [batch_correction_harmony.py](../panpipes/python_scripts/batch_correct_harmony.py)
This method is called in the [integration workflow](../panpipes/panpipes/pipeline_integration.py):

```
# rna HARMONY (see lines 146:184)
@follows(set_up_dirs)
@active_if(PARAMS['rna_run'])
@active_if(PARAMS['rna_tools'] is not None and 'harmony' in PARAMS['rna_tools'])
@originate("batch_correction/umap_rna_harmony.csv")
def run_harmony(outfile):
    cmd = """python %(py_path)s/batch_correct_harmony.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --integration_col %(rna_column)s
     --n_threads %(resources_threads_high)s
     --modality rna
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    harmony_params = PARAMS['rna']['harmony']
    if harmony_params['npcs'] is not None:
        cmd += " --harmony_npcs %s" % harmony_params['npcs']
    if harmony_params['sigma'] is not None:
        cmd += " --sigma_val %s" % harmony_params['sigma'] 
    if harmony_params['theta'] is not None:
        cmd += " --theta_val %s" % harmony_params['theta']       
    neighbor_params = PARAMS['rna']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/rna_harmony.log " 
     #job arguments
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)

```

The first 4 lines of the given code are ruffus decorators, which is what allows the definition of the tasks execution flow and conditions for the subsequent tasks in the data analysis pipeline. 

Here, we are specifying that `run_harmony` will run after the previous task `set_up_dirs` is complete run (see [line 31](../panpipes/panpipes/pipeline_integration.py#L31)). 
We're also specifying that the `run_harmony` task will run only if in the configuration file the `rna_run: ` voice is set to True and if `harmony` is one of the tools listed under `rna_tools: `, amongst those that can perform batch correction on rna that we include.
Finally the `@originate` decorator is what allows to generate the output of this task. 
That means that when running the `integration` workflow, the workflow will check if a file `"batch_correction/umap_rna_harmony.csv"` exists, if not, it will run the harmony integration task.

As you can see from the subsequent few lines, we are passing the parameters to the harmony python script which reads them leveragint the `argparse` library.
Let's look more closely at the [batch_correction_harmony.py](../panpipes/python_scripts/batch_correct_harmony.py) script to understand how `harmony` is run.
This scripts expects as input a MuData object with a layer called `rna`, on which PCA has been calculated (this happens in the `preprocess` workflow). 
After the necessary checks are carried out, the script calls `ho = hm.run_harmony(...)` which is the `harmony` implementation provided in the `harmonypy` package.

The new dimensionality reduction is computed and added in the relevant slot of the `rna` layer. Then, the knn graph is computed and a umap is also generated and added in the `.obsm` slot. Finally, the outputs are saved, namely a `MuData` with the new harmony representation and a csv file with the single cell umap coordinates we have computed.

While running this task, panpipes will print on screen a message:

```
2023-08-29 15:34:58,651 INFO main task - Task enters queue = 'pipeline_integration.run_harmony'
2023-08-29 15:34:59,228 INFO main execution - job-options: -J umap_rna_harmony.csv --constraint=skl --cpus-per-task=6 --partition=short
2023-08-29 15:34:59,228 INFO main execution - running statement: \
                                             python panpipes/python_scripts/batch_correct_harmony.py       --input_anndata pbmc_seuratv4.h5mu      --output_csv batch_correction/umap_rna_harmony.csv       --integration_col sample_id      --n_threads 6      --modality rna       --harmony_npcs 30 --sigma_val 0.1 --theta_val 1.0 --neighbors_method scanpy --neighbors_metric euclidean --neighbors_n_pcs 30 --neighbors_k 30 > logs/rna_harmony.log
2023-08-29 15:34:59,240 INFO main execution - job has been submitted with job_id 24942982

```


When the task is completed (so when the `"batch_correction/umap_rna_harmony.csv"` file is produced), the pipeline exits the task with a message

```
2023-08-29 15:41:09,811 INFO main task - Completed Task = 'pipeline_integration.run_harmony' 
```
and it's ready for the next task.

### 2. You can contribute too!

Now you know most of the important bits to understand how `panpipes` orchestrates the analysis flows. If you want to contribute a method to panpipes, you can make a fork of the github repo and open a Pull Request when your code is ready and tested. 
All you have to do is write a script like the one we showed before and make sure you provide its required parameters in the configuration file. 

You will need:

1. A script calling the necessary modules to run your analysis. Use a library like `argparse` to allow parsing the arguments to your script.
2. If the script you're calling needs a new parameter to be configured in the configuration .yml, remember to add it there, so the pipeline can read it when parsing the `PARAMS` argument.
   
*We recommend testing the script as a standalone before adding it as a formal component to the pipeline, so you can provide the inputs directly instead of waiting for the pipeline to run the previous tasks*

3. Add the script to the relevant workflow. For example, if you are adding a new multimodal integration method, the best place to add it is the `integration` workflow.
Check the decorators on the other tasks that are similar to the one you want to add, and make sure you're following the same flow principles where applicable. Remember that an analysis pipeline requires you to think in terms of inputs and outputs, and that tasks that accept the same input should be parallelized wherever possible.
For more guidance on ruffus decorators please check the [main guidelines](http://www.ruffus.org.uk/index.html) 


We are mostly python developers and we designed panpipes to be part of the scverse, so most of our code is written in python. However, if you would like to implement a component that uses R, check out the interoperability section in the [scverse documentation](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/scverse_data_interoperability.html).


The interaction with the Workflow Management System will be largely the same, but for an example of how we are already implementing R components please check the `clustree` component in the [`clustering` workflow](../panpipes/panpipes/pipeline_clustering.py#L223-L236).


We're happy to help, just open issues on our github page!









