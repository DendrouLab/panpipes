Normalization methods
=======================


panpipes currently supports the following normalization methods:

## RNA

1. Standard normalization using scanpy's [normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html) and [log1p](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html).

Additionally, the normalized data can be scaled using scanpy's [sc.pp.scale](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.scale.html)

## PROT

1. clr using [muon's prot processing](https://muon.readthedocs.io/en/latest/api/generated/muon.prot.pp.html), with the option to specify margin for normalization,
   1. clr margin= 0, normalize within each feature's distribution, across all cells
   2. clr margin= 1, normalize within each cells' counts distribution, across all features
   
    *if you come from R, please note that the [margins are transposed](https://en.wikipedia.org/wiki/The_Scream#/media/File:Edvard_Munch,_1893,_The_Scream,_oil,_tempera_and_pastel_on_cardboard,_91_x_73_cm,_National_Gallery_of_Norway.jpg) in the Python and anndata world*

2. dsb using [muon's prot processing](https://muon.readthedocs.io/en/latest/api/generated/muon.prot.pp.html)

## ATAC

1. Standard normalization using scanpy's [normalize_total]() and [log1p](). 
2. TFIDF with 3 flavours
   1. "signac", following [signac's defaults](https://stuartlab.org/signac/articles/pbmc_vignette#normalization-and-linear-dimensional-reduction).
    using [muon's atac processing](https://muon.readthedocs.io/en/latest/api/generated/muon.atac.pp.tfidf.html#muon.atac.pp.tfidf) 
   2. logTF: logging the TF term using using [muon's atac processing](https://muon.readthedocs.io/en/latest/api/generated/muon.atac.pp.tfidf.html#muon.atac.pp.tfidf) 
   3. logIDF: logging the IDF term using using [muon's atac processing](https://muon.readthedocs.io/en/latest/api/generated/muon.atac.pp.tfidf.html#muon.atac.pp.tfidf) 



Layers nomenclature
============================

Normalised data is saved in each modality slot in their specific layers:

    atac.layers["logTF_norm"] = atac.X.copy()
    

| method                 | layer           | modality      |
| ---------------------- | --------------- | ------------- |
| raw counts             | "raw_counts     | RNA/ATAC/PROT |
| standard log1p         | "logged_counts" | RNA or ATAC   |
| scaled counts          | "scaled_counts" | RNA or ATAC   |
| clr                    | "clr"           | PROT          |
| dsb                    | "dsb"           | PROT          |
| TFIDF (signac flavour) | "signac_norm"   | ATAC          |
| TFIDF (logTF)          | "logTF_norm"    | ATAC          |
| TFIDF (logIDF)         | "logIDF_norm"   | ATAC          |



