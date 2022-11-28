from anndata import AnnData
# from helpers import gen_adata
import panpipes.funcs as pnp
from scipy import sparse
import pytest
import pandas as pd
import numpy as np


@pytest.fixture()
def adata(nn=4):
    yield AnnData(sparse.eye(nn), 
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(nn)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(nn)]))


@pytest.fixture()
def mdata():
    yield MuData(
        {
            "mod1": AnnData(np.arange(0, 100, 0.1).reshape(-1, 10)),
            "mod2": AnnData(np.arange(101, 2101, 1).reshape(-1, 20)),
        }
    )

@pytest.fixture()
def anndata():
    yield AnnData(np.arange(50, 100, 1).reshape(-1, 10), 
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(5)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(10)]))

# def test_merge_with_adata_obs(adata):
#     adata['sample_id'] = 'a'
#     obs = pd.DataFrame(pd.DataFrame(index=[f"cell{i}" for i in range(4)]), data={'batch'==[0,0,1,1]})
