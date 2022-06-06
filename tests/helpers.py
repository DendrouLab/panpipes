from scipy import sparse
import numpy as np
import pandas as pd
from typing import Tuple
from anndata import AnnData

def gen_adata(
    shape: Tuple[int, int],
    X_type=sparse.csr_matrix,
    X_dtype=np.float32,
    # obs_dtypes,
    # var_dtypes,
    obsm_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
    varm_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
    layers_types: "Collection[Type]" = (sparse.csr_matrix, np.ndarray, pd.DataFrame),
) -> AnnData:
    """\
    Borrowed from https://github.com/theislab/anndata/blob/9520e5967ca0a12de92cd4672e185a796a9289e8/anndata/tests/helpers.py#L65
    Helper function to generate a random AnnData for testing purposes.
    Note: For `obsm_types`, `varm_types`, and `layers_types` these currently
    just filter already created objects.
    In future, these should choose which objects are created.
    Params
    ------
    shape
        What shape you want the anndata to be.
    X_type
        What kind of container should `X` be? This will be called on a randomly
        generated 2d array.
    X_dtype
        What should the dtype of the `.X` container be?
    obsm_types
        What kinds of containers should be in `.obsm`?
    varm_types
        What kinds of containers should be in `.varm`?
    layers_types
        What kinds of containers should be in `.layers`?
    """
    M, N = shape
    obs_names = pd.Index(f"cell{i}" for i in range(shape[0]))
    var_names = pd.Index(f"gene{i}" for i in range(shape[1]))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    # For #147
    obs.rename(columns=dict(cat="obs_cat"), inplace=True)
    var.rename(columns=dict(cat="var_cat"), inplace=True)

    if X_type is None:
        X = None
    else:
        X = X_type(np.random.binomial(100, 0.005, (M, N)).astype(X_dtype))
    obsm = dict(
        array=np.random.random((M, 50)),
        sparse=sparse.random(M, 100, format="csr"),
        df=gen_typed_df(M, obs_names),
    )
    obsm = {k: v for k, v in obsm.items() if type(v) in obsm_types}
    varm = dict(
        array=np.random.random((N, 50)),
        sparse=sparse.random(N, 100, format="csr"),
        df=gen_typed_df(N, var_names),
    )
    varm = {k: v for k, v in varm.items() if type(v) in varm_types}
    layers = dict(
        array=np.random.random((M, N)), sparse=sparse.random(M, N, format="csr")
    )
    layers = {k: v for k, v in layers.items() if type(v) in layers_types}
    obsp = dict(
        array=np.random.random((M, M)), sparse=sparse.random(M, M, format="csr")
    )
    varp = dict(
        array=np.random.random((N, N)), sparse=sparse.random(N, N, format="csr")
    )
    uns = dict(
        O_recarray=gen_vstr_recarray(N, 5),
        nested=dict(
            scalar_str="str",
            scalar_int=42,
            scalar_float=3.0,
            nested_further=dict(array=np.arange(5)),
        ),
        # U_recarray=gen_vstr_recarray(N, 5, "U4")
    )
    adata = AnnData(
        X=X,
        obs=obs,
        var=var,
        obsm=obsm,
        varm=varm,
        layers=layers,
        obsp=obsp,
        varp=varp,
        dtype=X_dtype,
        uns=uns,
    )
    return adata
