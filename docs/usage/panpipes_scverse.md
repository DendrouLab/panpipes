# Panpipes and the scverse


Panpipes is an analytical pipeline implemented in `python`, which orchestrates single cell analyses in an automated fashion. Panpipes relies on the `scverse`, an ecosystem of tools for single cell 'omics analysis in python.
Checkout the [scverse](https://scverse.org/) documentation for the packages included and how to contribute your package so we can support it in `panpipes` too! 



Panpipes has at its core [`AnnData`](https://anndata.readthedocs.io/en/latest/) and [`MuData`](https://mudata.readthedocs.io/en/latest/), for handling annotated data matrices in memory and on disk, with the best of the pandas and xarray functionalities.

<figure>
    <img src="https://github.com/DendrouLab/panpipes/blob/main/docs/img/anndata_schema.svg?raw=true" alt="img1" width="40%">
    <figcaption>`AnnData` is anndata is a container for handling annotated data matrices objects.</figcaption>
</figure>


<figure>
    <img src="https://github.com/DendrouLab/panpipes/blob/main/docs/img/mudata_paper.svg?raw=true" alt="img2" width="40%">
    <figcaption>`MuData` is a dictionary of `AnnData` objects.</figcaption>
</figure>

The workhorses for panpipes are `scanpy`, `muon` and `squidpy`, frameworks for analyzing single-cell gene expression, multimodal data and spatial trascriptomics.

For deep-learning based methods for uni or multimodal integration, we leverage the functionalities of `scvi-tools`, a library developing probabilistic models for single-cell omics data in PyTorch.

To get help with `panpipes`, you can open an issue on our [Github page](). 
Please use the [scverse discourse](https://discourse.scverse.org/) to document issues with `scverse` packages and get the help of other scverse users! 


Please use these links to familiarize with these data structures and frameworks:

- [AnnData](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [MuData Quickstart](https://muon.readthedocs.io/en/latest/notebooks/quickstart_mudata.html)
- [Single cell analysis with scanpy](https://scanpy.readthedocs.io/en/latest/) 
- [Multimodal analyses with muon](https://muon-tutorials.readthedocs.io/en/latest/)
- [Spatial analyses with squidpy](https://squidpy.readthedocs.io/en/stable/)
- [Scvi-tools](https://scvi-tools.org/)
- [Interoperability](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/scverse_data_interoperability.html)
- [Best-Practices for single cell multimodal analyses ](https://www.sc-best-practices.org/preamble.html)

