import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import muon as mu
from cgatcore import pipeline as P
import panpipes.funcs as pp
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice, X_is_raw

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# load arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--scaled_anndata',
                    default='adata_scaled.h5ad',
                    help='a preprocessed mudata/anndata')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_wnn.csv',
                    help='')
parser.add_argument('--figdir', default='./figures',
                    help='')
parser.add_argument('--n_neighbors',
                    help="number of neighbors", default=50)
parser.add_argument('--n_bandwidth_neighbors', default="",
                    help="")
parser.add_argument('--n_multineighbors', default=30,
                    help="neighbors k")
parser.add_argument('--metric',default="euclidean",
                    help="neighbor metric, e.g. euclidean or cosine")
parser.add_argument('--low_memory',default=True,
                    help="set to True by default if cells in dataset >50k")


args, opt = parser.parse_known_args()

L.info(args)
# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir

params = pp.io.read_yaml("pipeline.yml")
threads_available = multiprocessing.cpu_count()


if params['multimodal']['WNN']['modalities'] is not None:
    modalities= params['multimodal']['WNN']['modalities']
    modalities = [x.strip() for x in modalities.split(",")]
    L.info(f"using modalities :{modalities}")

L.info("running with batch corrections:")
wnn_params_bc = params['multimodal']['WNN']['batch_corrected'] 
if modalities is not None:
    wnn_params_bc= {k: wnn_params_bc[k] for k in wnn_params_bc.keys() & modalities}
L.info( wnn_params_bc )

mdata = mu.read(args.scaled_anndata)

if all(x in modalities for x in mdata.mod.keys()):
    tmp = mdata.copy()
    removed_mods = None
else:
    tmp = mdata.copy()
    removed_mods = list(set(mdata.mod.keys()) - set(modalities))
    L.info(f"removing modalities {removed_mods}")
    for rmod in removed_mods:
        del tmp.mod[rmod]

L.info("intersecting modality obs before running wnn")
mu.pp.intersect_obs(tmp)


#one could also check of obsp is not empty and use precomputed connectivities but i assume that if batch corrected the obsp is populated, if not, calc on the flight on pca
dict_graph = {}
for x in wnn_params_bc.keys():
    dict_graph[x] = {}
    if wnn_params_bc[x] is not None:
        dict_graph[x]["obsm"] = "X_" + wnn_params_bc[x]
        dict_graph[x]["anndata"] = "tmp/" + wnn_params_bc[x].lower() + "_scaled_adata_" + x +".h5ad"
    else: 
        dict_graph[x]["obsm"] = None

L.debug(dict_graph)

for kmod in dict_graph.keys():
    L.info(kmod)
    pkmod=params['multimodal']['WNN']['knn'][kmod]
    if dict_graph[kmod]["obsm"] is not None:
        if dict_graph[kmod]["obsm"] not in tmp.mod[kmod].obsm.keys():
            L.info("provided mdata doesn't have the desired obsm, just checking if it's bbknn you want.")
            if dict_graph[kmod]["obsm"] == "X_bbknn":
                if len(tmp.mod[kmod].obsp.keys()) > 0 and "neighbors" in tmp.mod[kmod].uns.keys() : 
                     L.info("i found a populated obsp slot and I assume it's bbknn")
                else:
                    if dict_graph[kmod]["anndata"] is not None:
                        L.info("reading precomputed connectivities for bbknn")
                        adata = mu.read(dict_graph[kmod]["anndata"])
                        tmp.mod[kmod].obsp = adata.obsp.copy()
                        tmp.mod[kmod].obsm[dict_graph[kmod]["obsm"]] = adata.obsm[dict_graph[kmod]["obsm"]].copy()
                        tmp.mod[kmod].uns["neighbors"]= adata.uns["neighbors"].copy()
                        tmp.update()
            else:
                if dict_graph[kmod]["anndata"] is not None:
                    L.info("provided mdata doesn't have the desired obsm. reading the batch corrected data from another stored object")
                    adata = mu.read(dict_graph[kmod]["anndata"])
                    L.debug(kmod + "object")
                    L.debug(adata)
                    tmp.mod[kmod].obsm[dict_graph[kmod]["obsm"]] = adata.obsm[dict_graph[kmod]["obsm"]].copy()
                    tmp.mod[kmod].obsp = adata.obsp.copy()
                    tmp.mod[kmod].uns['neighbors'] = adata.uns['neighbors'].copy()
                    tmp.update()
                    repuse = dict_graph[kmod]["obsm"] 
                else:
                    L.info("could not find the desired obsm and the anndata slot is empty, will calculate on the flight")
                    if kmod =="atac":
                        if "X_lsi" in tmp.mod[kmod].obsm.keys():
                            repuse = "X_lsi"
                        else:
                            repuse = "X_pca"
                        L.info("falling back on %s" %(repuse) )
            
                    L.info("calculating neighbours")
                    if repuse != "X_bbknn":
                        run_neighbors_method_choice(tmp.mod[kmod], 
                            method=pkmod['method'], 
                            n_neighbors=int(pkmod['k']), 
                            n_pcs=min(int(pkmod['npcs']), mdata.var.shape[0]-1), #this should be the # rows of var, not obs
                            metric=pkmod['metric'], 
                            #does this throw an error if no PCA for any single mod is stored?
                            use_rep=repuse,
                            nthreads=max([threads_available, 6]))
        else:
            L.info("Using %s" %(dict_graph[kmod]["obsm"]))            
    else:
        L.info("could not find the desired obsm and the anndata slot is empty, will calculate on the flight")
        repuse ="X_pca"
        if kmod =="atac":
            if "X_lsi" in tmp.mod[kmod].obsm.keys():
                repuse = "X_lsi"
            else:
                repuse = "X_pca"
            L.info("falling back on %s" %(repuse) )

        L.info("calculating neighbours")
        
        run_neighbors_method_choice(tmp.mod[kmod], 
            method=pkmod['method'], 
            n_neighbors=int(pkmod['k']), 
            n_pcs=min(int(pkmod['npcs']), mdata.var.shape[0]-1), #this should be the # rows of var, not obs
            metric=pkmod['metric'], 
            #does this throw an error if no PCA for any single mod is stored?
            use_rep=repuse,
            nthreads=max([threads_available, 6]))

L.debug(tmp)
tmp.update()
L.debug(tmp)
L.info("Now running WNN")

mu.pp.neighbors(tmp, 
                n_neighbors= int(args.n_neighbors),
                n_bandwidth_neighbors= int(args.n_bandwidth_neighbors),
                n_multineighbors= int(args.n_multineighbors),
                metric= args.metric,
                low_memory= check_for_bool(args.low_memory),   
                key_added='wnn')
L.info("WNN finished now calculate umap")

mu.tl.umap(tmp,min_dist=0.4, neighbors_key='wnn')
#For taking use of mdata.obsp['connectivities'], it’s scanpy.tl.leiden() that should be used. not muon.tl.leiden
sc.tl.leiden(tmp, neighbors_key='wnn', key_added='leiden_wnn') 

umap = pd.DataFrame(tmp.obsm['X_umap'], tmp.obs.index)
umap.to_csv(args.output_csv)

#add back modalities that were not used for this run

if removed_mods is not None:
    for rmd in removed_mods:
        tmp.mod[rmd] = mdata.mod[rmd].copy()

tmp.write("tmp/wnn_scaled_adata.h5mu")

L.info("Done")

