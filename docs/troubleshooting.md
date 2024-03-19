# Troubleshooting common issues when running panpipes

### What to do when the pipeline breaks mid-run

Sometimes the pipeline will stop because, for example,  a parameter is wrong in the config file, or a path is not
accurate in the YAML file. 

TODO add examples of directories when pipeline breaks mid-run

**Solution**
First: check the log files to see what went wrong.
- Inspect the pipeline.log file: the bottom of the file will print the error that broke the pipeline.
TODO: add text from a failed pipeline.log as an example
- In this case the pipeline failed at ...: inspect the log file for this process in logs/xxx.log

You can fix the issue by ...

Second: Before re-running panpipes, we recommend deleting any intermediate files that were created in the previous run which
broke halfway through.  



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

