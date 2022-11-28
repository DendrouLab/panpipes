import pytest
# from scpipelines.funcs.processing import extract_parameter_from_fname
import panpipes.funcs as pnp
from anndata import AnnData
from muon import MuData
import pandas as pd
import numpy as np

@pytest.fixture()
def mudata():
    yield MuData(
        {
            "mod1": AnnData(np.arange(0, 100, 0.1).reshape(-1, 10)),
            "mod2": AnnData(np.arange(101, 2101, 1).reshape(-1, 20)),
        }
    )

@pytest.fixture()
def anndata():
    yield AnnData(np.arange(0, 50, 1).reshape(-1, 10), 
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(5)],
                         data={'sample_id': ['a', 'a', 'b', 'b', 'c']}),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(10)]))



@pytest.fixture()
def anndata_with_obs():
    yield AnnData(np.arange(0, 50, 1).reshape(-1, 10), 
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(5)],
                         data={'sample_id': ['a', 'a', 'b', 'b', 'c'],
                               'batch': [1,1,2,2,3]}),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(10)])
    )

@pytest.fixture()
def anndata_with_var():
    yield AnnData(np.arange(0, 50, 1).reshape(-1, 10), 
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(5)],
                         data={'sample_id': ['a', 'a', 'b', 'b', 'c']}),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(10)],
                         data={'feature_type': "Gene Expression",
                              'new_index': [f"newgene{i}" for i in range(10)]}
                        )
    )


def test_extract_parameter_from_fname():
    assert pnp.pp.extract_parameter_from_fname("test_res0.6_cluster.txt.gz", "res", "test") == 0.6
    assert pnp.pp.extract_parameter_from_fname("test_res1.0_cluster.txt.gz", "res", "test") == 1
    assert pnp.pp.extract_parameter_from_fname("test_methodeuclidean_cluster.txt.gz", "method", "test") == "euclidean"
    assert pnp.pp.extract_parameter_from_fname("method_methodeuclidean_cluster.txt.gz", "method", "method") == "euclidean"

def test_is_float_try():
    assert pnp.pp.is_float_try("1") is True
    assert pnp.pp.is_float_try("cheese") is False

def test_splitall():
    assert pnp.pp.splitall("path") == ['path']
    assert pnp.pp.splitall("path/to/cheese") == ['path', 'to', 'cheese']
    assert pnp.pp.splitall("/path/to/cheese") == ['path', 'to', 'cheese']

def test_test_file_or_value():
    assert pnp.pp.test_file_or_value(4) == "value"
    assert pnp.pp.test_file_or_value(__file__) == "file"
    with pytest.raises(ValueError):
        pnp.pp.test_file_or_value("cheese")

def test_which_ind():
    assert pnp.pp.which_ind([True, True, False]) == [0,1]
    assert pnp.pp.which_val([True, True, False], [1,2,3]) == [1,2]

def test_check_for_bool():
    assert pnp.pp.check_for_bool('True') is True
    assert pnp.pp.check_for_bool(True) is True
    assert pnp.pp.check_for_bool('False') is False
    assert pnp.pp.check_for_bool(False) is False
    with pytest.raises(TypeError):
        pnp.pp.check_for_bool('cheese')
    with pytest.raises(TypeError):
        pnp.pp.check_for_bool(4)

def test_intersection():
    assert pnp.pp.intersection(['a', 'b', 'c'], ['c', 'd', 'e']) == ['c']


def assert_match_anndata(ad1, ad2):
    assert np.array_equal(ad1.X, ad2.X)
    assert all(ad1.obs == ad2.obs)
    assert all(ad1.var == ad2.var)

def test_concat_adatas(anndata):
    ad1 = anndata
    ad1 = ad1.copy()
    ad1.obs['sample_id'] = 'a'
    ad2 = anndata
    ad2 = ad2.copy()
    ad2.obs['sample_id'] = 'b'
    # this is what the merged version should look like :-) 
    adata_full = AnnData(
        np.concatenate([np.arange(0, 50, 1),np.arange(0, 50, 1)]).reshape(-1, 10),
                obs = pd.DataFrame(index=['cell0-a', 'cell1-a', 'cell2-a', 'cell3-a', 'cell4-a', 
                                            'cell0-b', 'cell1-b','cell2-b', 'cell3-b', 'cell4-b'],
                                    data={'sample_id': ['a', 'a', 'a', 'a', 'a', 
                                                        'b', 'b', 'b', 'b', 'b']}, dtype="category"),
                var=pd.DataFrame(index=[f"gene{i}" for i in range(10)])
    )
    merg_data = pnp.pp.concat_adatas([ad1, ad2],
        batch_key="sample_id", 
        batch_categories=['a', 'b'], 
        join_type="inner")
    # testing the case where there are two anndatas to be merged
    assert np.array_equal(merg_data.X, adata_full.X)
    assert_match_anndata(merg_data, adata_full)
    # testing the case where there is only 1 anndata sent to function
    assert_match_anndata(pnp.pp.concat_adatas([ad1], batch_key="sample_id", 
        batch_categories=['a']), ad1)
    # test the corner case where an anddata is passed
    with pytest.raises(TypeError):
        pnp.pp.concat_adatas(ad1, batch_key="sample_id", 
        batch_categories=['a'])

def test_merge_with_adata_obs(anndata, anndata_with_obs):
    new_obs = pd.DataFrame(data={'sample_id': ['a', 'b', 'c'],
                             'batch': [1, 2, 3]})
    # test exceptions
    with pytest.raises(TypeError):
        pnp.pp.merge_with_adata_obs("cheese", new_obs, "sample_id")
    with pytest.raises(TypeError):
        pnp.pp.merge_with_adata_obs(anndata, "cheese", "sample_id")
    with pytest.raises(KeyError):
        pnp.pp.merge_with_adata_obs(anndata, new_obs, "cheese")
    with pytest.raises(KeyError):
        pnp.pp.merge_with_adata_obs(anndata, new_obs, "batch")
    anndata.obs = pnp.pp.merge_with_adata_obs(anndata, new_obs, "sample_id")
    assert_match_anndata(anndata, anndata_with_obs)


def test_merge_with_adata_obs_inplace(anndata, anndata_with_obs):
    new_obs = pd.DataFrame(data={'sample_id': ['a', 'b', 'c'],
                             'batch': [1, 2, 3]})
    pnp.pp.merge_with_adata_obs(anndata, new_obs, "sample_id", inplace=True)
    assert_match_anndata(anndata, anndata_with_obs)

def test_remove_unused_categories():
    df = pd.DataFrame(data={"a1":['a', 'b', 'c'], 'a2':['c', 'd', 'e']}, dtype='category')
    df = df.iloc[0:2,:]
    pnp.pp.remove_unused_categories(df)
    assert df['a1'].cat.categories.tolist() == ['a', 'b']
    assert df['a2'].cat.categories.tolist() == ['c', 'd']
    with pytest.raises(TypeError):
        pnp.pp.remove_unused_categories("cheese")

def test_add_var_mtd(anndata, anndata_with_var):
    new_var = pd.DataFrame(index=[f"gene{i}" for i in range(10)],
                         data={'feature_type': "Gene Expression",
                                'new_index': [f"newgene{i}" for i in range(10)]}
                        )
    pnp.pp.add_var_mtd(anndata, new_var, left_on="index", right_on="new_index")
    assert_match_anndata(anndata, anndata_with_var)


def test_add_var_mtd_types(anndata):
    with pytest.raises(TypeError):
        pnp.pp.add_var_mtd(anndata, 'cheese')
    new_var = pd.DataFrame(index=[f"gene{i}" for i in range(10)],
                         data={'feature_type': "Gene Expression",
                                'new_index': [f"newgene{i}" for i in range(10)]}
                        )
    with pytest.raises(TypeError):
        pnp.pp.add_var_mtd('cheese', new_var)
    # check mismatch var df raises coorect errors


# def test_add_var_mtd_warnings(anndata, anndata_with_var):
#     ad = anndata.copy()
#     new_var = pd.DataFrame(index=[f"gene{i}" for i in range(9)],
#                          data={'feature_type': "Gene Expression",
#                                 'new_index': [f"newgene{i}" for i in range(9)]})
#     # with pytest.warns(UserWarning):
#     #     pnp.pp.add_var_mtd(ad, new_var)


def test_update_var_index(anndata):
    new_index = [f"newgene{i}" for i in range(10)]
    anndata.var['new_index'] = new_index
    pnp.pp.update_var_index(anndata, "new_index")
    assert all(anndata.var.index == new_index)


