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

