<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>
# Integration methods implemented in panpipes

The panpipes integration pipeline implements a variety of tools to batch correct individual modalities and/or integrate across modalities to produce a reduced dimension representation of the dataset.<br>
There are different tools available for each modality such as RNA (also referred to as GEX), PROT (can be referred to as ADT) and ATAC which can be run as required before running `panpipes integration make merge_batch_correction`
ro create the final object with the reduced dimension represented.<br> 

The ideal way to run `panpipes integration` is to use the output `MuData`file from `panpipes preprocess` since it will already be in the required format. 
However, if using independent MuData the object should contain normalised data in the X slot of each modality, a ‘raw_counts’ layer in each modality, and a sample_id column in each slot of the obs and the outer obs.

| method       | type of integration | modalities | code 
| ------------ | ------------------- | ------------- |
| raw counts   | "raw_counts"        | RNA/ATAC/PROT |
