Reference Mapping (refmap)
======

The reference mapping pipeline `panpipes refmap` implements scvi-tools reference mapping tools `scvi`, `totalvi` and `scanvi`, and the Muon implementation of `wnn`.
For the scvi-tools methods, you can supply an anndata/mudata containing RNA data, or the simply the saved model data from a previous run. You can even use an scvi model created previously by `panpipes integration`.

Steps:
------

1.  In a new folder, generate config file for integration,
    `panpipes refmap config` and edit the pipeline.yml file.
2.  Run complete refmap pipeline with `panpipes refmap make full`


## Expected structure of MuData object
The ideal way to run `panpipes refmap` is to use the output mudata file from `panpipes preprocess`, as this will make sure the MuData object has correctly names layers and slots. 