#Troubleshooting common issues when running panpipes

### What to do when the pipeline breaks mid run

Sometimes the pipeline will stop because, for example,  a parameter is wrong in the config file, or a path was not
accurately changed in the yaml file. 

First: check the log file to see what went wrong and fix the issue 

Second: Before re running make sure to delete any intermediate files that where created in the previous run which
broke halfway through, to ensure that you can fully reattempt.  
