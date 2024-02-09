#Troubleshooting common issues when running panpipes


### What to do when the pipeline breaks mid run

Sometimes the pipeline will stop because, for example,  a parameter is wrong in the config file, or a path was not
accurately changed in the yml file. 

Firslty: check the log file to see what went wrong and fix the issue 

Secondly: Before re running make sure to delete any intermediate files that where created in the previous run which
broke halfway throught, to ensure that you can fully reattempt.  
