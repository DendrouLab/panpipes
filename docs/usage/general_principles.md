

General principles for running pipelines
----------------------------------------

Run the pipeline from the login node on your server, it will use in
built the job submission system to submit jobs. To prevent jobs being cancelled when you log out, we recommend using either the `nohup` command, or a terminal multiplexer such a [tmux](https://github.com/tmux/tmux/wiki) or [screen](https://www.gnu.org/software/screen/manual/screen.html).

Navigate to the directory where you want to run your analysis (this
should not be within the panpipes folder, or your virutal environment
folder)

    mkdir data_dir/
    cd data_dir/
    panpipes qc_mm config

This will produce two files, `pipeline.log` and `pipeline.yml`

Edit `pipeline.yml` as appropriate for your data, following the
instructions within the yml file.

Then check which jobs will run with the command

    panpipes qc_mm show full

The output of this will show a list of tasks that will be run as part of
the pipeline.

To run use the command

    panpipes qc_mm make full

Occasionally you might want to run tasks individually (e.g. to assess
outputs before deciding the parameters for the next step) In order to do
this you can run any task in the `show full` list such as:

    panpipes qc_mm make plot_tenx_metrics



Creating a [sample submission file](setup_for_qc_mm) for the QC pipeline.
