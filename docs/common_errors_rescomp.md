
# "pipeline exits with 0 errors":
This means something has gone wrong with one of the sub-scripts, rather than the main pipeline scripts. 
Try to work out which step it failed on, and check the approapriate log files. 
#### Solutions:
- 45% of the time this error occurs when one of the inputs or the pipeline.yml is wrongly formatted.
- 45% of the time this error occurs when one of the required libraries is not installed.
- 5% of the time code is actually broken, because you've found a corner case that hasn't been accounted for in the code.
- Final 5% magically works when resubmitted, ¯\\_(ツ)_/¯

# "Illegal Instruction" errors:
You are trying to run code using software that has been compiled for a different node.
Occurs when code compiuled on rescomp1 (skylake) gets submitted to ivybridge.
Solution: edit the .cgat.yml file to restrict which nodes the jobs get submitted to.

# Find markers jobs are being submitted sequentially rather than in parallel:
This is an inconsistent bug and I don't know how to fix it. 
Solution: If it occurs, just stop and restart the pipeline.


CRG 
2021-11-08
