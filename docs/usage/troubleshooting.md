# Troubleshooting common issues when running panpipes

### What to do when the pipeline breaks mid-run

Sometimes the pipeline will stop because, for example,  a parameter is wrong in the config file, or a path is not accurate in the YAML file. 

Let's look at an example:


```
2024-04-24 17:48:48,048 INFO running statement: \
#                              python /Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/python_scripts/batch_correct_none.py       --input_anndata teaseq.h5mu      --output_csv batch_correction/umap_rna_none.csv      --integration_col dataset       --neighbors_method scanpy --neighbors_metric euclidean --neighbors_n_pcs 30 --neighbors_k 30 > logs/rna_no_correct.log
# 2024-04-24 17:50:14,045 INFO {"task": "'pipeline_integration.run_no_batch_umap'", "task_status": "update", "task_total": 1, "task_completed": 0, "task_completed_percent": 0.0}
# 2024-04-24 17:50:14,045 ERROR  \
#                                   Exception #2 \
#                                     'builtins.OSError(--------------------------------------- \
#                                   Child was terminated by signal -1:  \
#                                   The stderr was:  \
#                                   Traceback (most recent call last): \
#                                     File "/Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/python_scripts/batch_correct_none.py", line 61, in <module> \
#                                       adata = mu.read(args.input_anndata +"/" + args.modality) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/mudata/_core/io.py", line 602, in read \
#                                       return read_h5ad(filepath, m[2], **kwargs) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/mudata/_core/io.py", line 563, in read_h5ad \
#                                       with h5py.File(filename, hdf5_mode) as f_root: \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/h5py/_hl/files.py", line 562, in __init__ \
#                                       fid = make_fid(name, mode, userblock_size, fapl, fcpl, swmr=swmr) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/h5py/_hl/files.py", line 235, in make_fid \
#                                       fid = h5f.open(name, flags, fapl=fapl) \
#                                     File "h5py/_objects.pyx", line 54, in h5py._objects.with_phil.wrapper \
#                                     File "h5py/_objects.pyx", line 55, in h5py._objects.with_phil.wrapper \
#                                     File "h5py/h5f.pyx", line 102, in h5py.h5f.open \
#                                   FileNotFoundError: [Errno 2] Unable to open file (unable to open file: name = 'teaseq.h5mu', errno = 2, error message = 'No such file or directory', flags = 0, o_flags = 0) \
#                                    \
#                                   python /Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/python_scripts/batch_correct_none.py       --input_anndata teaseq.h5mu      --output_csv batch_correction/umap_rna_none.csv      --integration_col dataset       --neighbors_method scanpy --neighbors_metric euclidean --neighbors_n_pcs 30 --neighbors_k 30 > logs/rna_no_correct.log \
#                                   -----------------------------------------)' raised in ... \
#                                      Task = def pipeline_integration.run_no_batch_umap(...): \
#                                      Job  = [None -> batch_correction/umap_rna_none.csv] \
#                                    \
#                                   Traceback (most recent call last): \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/ruffus/task.py", line 712, in run_pooled_job_without_exceptions \
#                                       return_value = job_wrapper(params, user_defined_work_func, \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/ruffus/task.py", line 608, in job_wrapper_output_files \
#                                       job_wrapper_io_files(params, user_defined_work_func, register_cleanup, touch_files_only, \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/ruffus/task.py", line 540, in job_wrapper_io_files \
#                                       ret_val = user_defined_work_func(*(params[1:])) \
#                                     File "/Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/panpipes/pipeline_integration.py", line 89, in run_no_batch_umap \
#                                       P.run(cmd, **job_kwargs) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/cgatcore/pipeline/execution.py", line 1244, in run \
#                                       benchmark_data = r.run(statement_list) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/cgatcore/pipeline/execution.py", line 1029, in run \
#                                       raise OSError( \
#                                   OSError: --------------------------------------- \
#                                   Child was terminated by signal -1:  \
#                                   The stderr was:  \
#                                   Traceback (most recent call last): \
#                                     File "/Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/python_scripts/batch_correct_none.py", line 61, in <module> \
#                                       adata = mu.read(args.input_anndata +"/" + args.modality) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/mudata/_core/io.py", line 602, in read \
#                                       return read_h5ad(filepath, m[2], **kwargs) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/mudata/_core/io.py", line 563, in read_h5ad \
#                                       with h5py.File(filename, hdf5_mode) as f_root: \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/h5py/_hl/files.py", line 562, in __init__ \
#                                       fid = make_fid(name, mode, userblock_size, fapl, fcpl, swmr=swmr) \
#                                     File "/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/lib/python3.10/site-packages/h5py/_hl/files.py", line 235, in make_fid \
#                                       fid = h5f.open(name, flags, fapl=fapl) \
#                                     File "h5py/_objects.pyx", line 54, in h5py._objects.with_phil.wrapper \
#                                     File "h5py/_objects.pyx", line 55, in h5py._objects.with_phil.wrapper \
#                                     File "h5py/h5f.pyx", line 102, in h5py.h5f.open \
#                                   FileNotFoundError: [Errno 2] Unable to open file (unable to open file: name = 'teaseq.h5mu', errno = 2, error message = 'No such file or directory', flags = 0, o_flags = 0) \
#                                    \
#                                   python /Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/python_scripts/batch_correct_none.py       --input_anndata teaseq.h5mu      --output_csv batch_correction/umap_rna_none.csv      --integration_col dataset       --neighbors_method scanpy --neighbors_metric euclidean --neighbors_n_pcs 30 --neighbors_k 30 > logs/rna_no_correct.log \
#                                   ----------------------------------------- \
#                                    \

```

**Solution**

- Inspect the pipeline.log file: the bottom of the file will print the error that broke the pipeline. In this case the error tells us that the file we expected to use for the integration, 
`teaseq.h5mu` is not in the directory


```
FileNotFoundError: [Errno 2] Unable to open file (unable to open file: name = 'teaseq.h5mu', errno = 2, error message = 'No such file or directory', flags = 0, o_flags = 0) \
#                                    \
```
Indeed, in the configuration file, we specified the wrong path to the file :

```
preprocessed_obj: ../preprocess/teaseq.h5mu
```
In this case, the issue can be fixed by specifying the right path.

Checking the main directory we can see that the pipeline has left its footprint when failing:


```
-rw-r--r--@  1 fx  2125895594   14405 Apr 22 17:40 teaseq.h5mu
-rw-r--r--@  1 fx  2125895594   14405 Apr 22 17:46 pipeline.yml
drwxr-xr-x   8 fx  2125895594     256 Apr 22 17:47 figures
drwxr-xr-x  20 fx  2125895594     640 Apr 22 17:47 logs
-rwxrwx---   1 fx  2125895594    1024 Apr 22 17:48 ctmp4d4q4k3l.sh
drwxr-xr-x   2 fx  2125895594      64 Apr 22 17:48 batch_correction
-rw-r--r--   1 fx  2125895594     410 Apr 22 17:50 ctmp4d4q4k3l.sh.times
-rw-r--r--   1 fx  2125895594  325715 Apr 22 17:50 pipeline.log

```

We can also inspect the bash script, to see what was the command that failed last. 

If multiple parallel processes were running and failed because all used the same input parameters, you would see multiple bash scripts.

```
#!/bin/bash -eu
set -o pipefail

cd teaseq

mkdir -p batch_correction
umask 002
mkdir -p /var/folders/q1/jqpy4lk14cg0950t1mfrcf7cw66yvr/T/ctmpvtnf0426
export TMPDIR=/var/folders/q1/jqpy4lk14cg0950t1mfrcf7cw66yvr/T/ctmpvtnf0426

clean_temp() { rm -rf /var/folders/q1/jqpy4lk14cg0950t1mfrcf7cw66yvr/T/ctmpvtnf0426; }

info() { echo 'benchmark'; hostname; times; }

clean_all() { clean_temp; info; }

trap clean_all EXIT

umask 002
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PATH=/Users/fabiola.curion/Documents/devel/miniconda3/envs/pipeline_env/bin:$PATH
python /Users/fabiola.curion/Documents/devel/github/panpipes/panpipes/python_scripts/batch_correct_none.py       --input_anndata teaseq.h5mu      --output_csv batch_correction/umap_rna_none.csv      --integration_col dataset       --neighbors_method scanpy --neighbors_metric euclidean --neighbors_n_pcs 30 --neighbors_k 30 > logs/rna_no_correct.log

```

Inspecting the log file shows that indeed the input file is missing.

```
2024-04-22 17:50:12,756: INFO - reading data and starting integration pipeline with script: 
2024-04-22 17:50:12,756: INFO - batch_correct_none.py
2024-04-22 17:20:13,710: INFO - Namespace(iput_anndata='teaseq.h5mu', output_csv='batch_correction/umap_rna_scvi.csv', integration_col='dataset', neighbors_method='scanpy', neighbors_metric='euclidean', neighbors_n_pcs='30',neighbors_k='30')
2024-04-22 17:20:23,710: INFO - missing input anndata
```



**Solution**: Fix the path in the pipeline.yml. Before re-running panpipes, we recommend deleting any intermediate files that were created in the previous run which broke halfway through.  



### Error in Plot 10x metric task in ingestion when not starting from CellRanger outputs
**NoneType error**

First: check the log files to see what went wrong.
- In this case the pipeline failed at:

```
TypeError: argument of type 'NoneType' is not iterable \
```
- After checking the error, check which log file to inspect by checking the Job line:
```
Task = def pipeline_ingest.aggregate_tenx_metrics_multi(...): \
Job  = [None -> logs/tenx_metrics_multi_aggregate.log] \
```

- Inspect the log file for this process in `logs/tenx_metrics_multi_aggregate.log`
- Check which Python script the code failed to know which task failed by looking at. In this case, the error was in:

  ```
  File "../panpipes/panpipes/python_scripts/aggregate_cellranger_summary_metrics.py"
  ```
- This indicates that in the pipeline.yml some parameter involving `10x metric` was set wrong, and by going to the beginning of the pipeline.log file it is clear that the `plot_10X_metrics` was set to `True` although the data did not come from CellRanger, which means `plot_10X_metrics` should be set to `False` since the input data was an Andata object.
- Thus, by changing the `pipeline.yml` `plot_10X_metrics` to `False` the issue is fixed
- To rerun, first delete the `logs/tenx_metrics_multi_aggregate.log` file
- Then run `panpipes ingest make full`



### Error in Plot 10x metric task in ingestion when starting from CellRanger outputs using 10X.h5 format
**KeyError**

First: check the log files to see what went wrong.
- In this case the pipeline failed at:

```
raise KeyError(key) from err \
KeyError: 'Median UMI counts per cell' \
```

- After checking the error, check what was the last task run:

```
Traceback (most recent call last): \
File "/data/leuven/344/vsc34406/Miniconda3/envs/Panpipes/lib/python3.9/site-packages/panpipes/python_scripts/aggregate_cellranger_summary_metrics.py", line 268, in <module> \
```

- In this case, the error was while running the Python script `aggregate_cellranger_summary_metrics.py`
- This means that a particular string could not be found. In this case, it was the `Median UMI counts per cell`






### Wrong cell_cycle gene list inputted
**No valid genes were passed for scoring**

First: check the log files to see what went wrong.
- In this case the pipeline failed at:
```
ValueError: No valid genes were passed for scoring. \
```

 - After checking the error, check which log file to inspect by checking the Job line:
```
 Task = def pipeline_ingest.run_rna_qc(...): \
 Job  = [None -> logs/run_scanpy_qc_rna.log, mouse_integrated_cell_metadata.tsv, mouse_integrated_unfilt.h5mu] \
```
- In this case, also check the WARNING:

```
  WARNING: genes are not in var_names and ignored: ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8'] \
```
- This shows that the gene list provided is not matching with the input data since **none** of the genes can be found. For example, this gene list is for humans only and the input data was from a mouse model which requires a mouse-specific gene list
- To fix the issue simply change the preprocess pipeline.yml `custom_genes_file` to a mouse gene list provided in the `panpipes/resources`
- Delete the log folders and any temp files
- Run `panpipes ingest make ful`


### Running preprocessed with RNA-only PlotQC error 
**No such file or directory**
First: check the log files to see what went wrong.
- In this case the pipeline failed at:
```
FileNotFoundError: [Errno 2] No such file or directory:
```
 - After checking the error, check which log file to inspect by checking the Job line:
```
 Task = def pipeline_preprocess.filter_mudata(...): \
 Job  = [None -> humanised_preprocessed.h5mu] \
```
- Before proceeding, double-check the `pipeline.yml` to see if the correct path was provided for the mudata object.
- Once the user is sure the path is correct, the next step is to check the modalities in the `pipeline.yml`
- In this case, the mudata object contains only RNA and that was the only modality set for `True`, however in the `plotqc` variable there are metrics for `prot_metrics` when there is no protein in the data.
- To fix the issue, clear the metrics in `prot_metrics`
- Delete the log folder and any temp files
- Re-run `panpipes preprocess make full`
