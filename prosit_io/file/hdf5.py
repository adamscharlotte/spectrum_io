import pandas as pd
import numpy as np
import h5py
from scipy.sparse import coo_matrix
from typing import Union, Optional, List
import scipy
import pandas as pd
import logging
import threading

logger = logging.getLogger(__name__)

META_DATA_KEY = 'meta_data'
INTENSITY_RAW_KEY = 'raw_intensity'
INTENSITY_PRED_KEY = 'pred_intensity'
MZ_RAW_KEY = 'raw_mz'

def read_file(path: str, key: str) -> pd.DataFrame:
    """
    Read hdf5 file and return dataframe with contents.
    With possibility to partial load for memory issues.
    :param path: The path to the hdf5 file to read.
    :param key: The key of the dataset/group of interest.
    :return: a pandas DataFrame
    """
    try:
        if key.startswith("sparse"):
            with h5py.File(path, 'r') as f:
                print(f.keys())
                values = f[f"{key}/values"]
                i = f[f"{key}/i"]
                j = f[f"{key}/j"]
                shape = f[f"{key}/shape"]
                sparse_data = coo_matrix((values, (i, j)), shape)
                df = pd.DataFrame.sparse.from_spmatrix(sparse_data)
                if f"{key}/column_names" in f.keys():
                    df.columns = f[f"{key}/column_names"]
                if f"{key}/index" in f.keys():
                    df.index = f[f"{key}/index"]
        else:
            df = pd.read_hdf(path, key=key)
        return df
    except Exception as e:
        logger.exception(e)

def thread_this(fn):
    def run(*args, **kwargs):
        t = threading.Thread(target=fn, args=args, kwargs=kwargs)
        t.start()
        return t
    return run

@thread_this
def write_file(
        data_sets:[Union[pd.DataFrame, scipy.sparse.spmatrix]],
        path: str,
        dataset_names: [str],
        column_names: [[List[str]]] = None,
        ):
    """
    Writes several datasets (spectra) to hdf5 file.
    :param data: list of datasets
    :param path: The path to store the file to.
    :param dataset_names: List of dataset names
    :param column_names: List of column_names
    :return:
    """
    index = 0
    for data_set,dataset_name in zip(data_sets,dataset_names):
        if isinstance(data_set, pd.DataFrame):
            write_dataset(data_set, path, dataset_name)
        else:
            write_dataset(data_set, path, dataset_name, mode='a', column_names= column_names[index])
            index += 1



def write_dataset(
    data:Union[pd.DataFrame, scipy.sparse.spmatrix],
    path: str,
    dataset_name: str,
    mode: str = 'w',
    compression: Optional[Union[str, bool]] = True,
    column_names: Optional[Union[List[str], None]] = None,
    index: Optional[Union[List[str], None]] = None
):
    """
    Writes or appends dataset to an hdf5 file.
    :params
        data: The data to store. Can be a pandas DataFrame or a scipy Sparsematrix.
        path: The path to store the file to.
        dataset_name: The key in the hdf5 file under which to store the data.
        mode: The method when writing the data. Use 'a' to append to an existing file or 'w' to overwrite.
        compression: Optional, the method for compressing data. Check pandas.DataFrame.to_hdf docs for supported
            compression methods in case of providing a pandas DataFrame and h5py.Dataset docs in case of providing
            a sparse matrix. If providing False, no compression is applied; if providing True, defaults to 'zlib'
            in case of providing a pandas DataFrame and 'gzip' if providing a sparse matrix.
            standard compression depending on the data type given. Default: True.
        column_names: Optional, additional column column_names. Ignored if providing a pandas DataFrame. Default: None.
        index: Optional, additional index. Ignored if providing a pandas DataFrame. Default: None.
    """
    if isinstance(compression, bool):
        if compression:
            if isinstance(compression, pd.DataFrame):
                compression = 'zlib'
            elif isinstance(compression, scipy.sparse.spmatrix):
                compression = 'gzip'
            else:
                compression = None
        else:
            compression = None
    try:
        if isinstance(data, pd.DataFrame):
            data.to_hdf(path, key=dataset_name, mode=mode, complib=compression)
        elif isinstance(data, scipy.sparse.spmatrix):
            with h5py.File(path, mode) as f:
                group_name = f"sparse_{dataset_name}"
                f.create_group(group_name)
                i,j, values = scipy.sparse.find(data)
                shape = data.shape
                f.create_dataset(f"{group_name}/i", data=i, compression=compression, dtype=int)
                f.create_dataset(f"{group_name}/j", data=j, compression=compression, dtype=int)
                f.create_dataset(f"{group_name}/values", data=values, compression=compression, dtype=float)
                f.create_dataset(f"{group_name}/shape", data=shape, shape=(2,), dtype=int)
                if column_names:
                    f.create_dataset(f"{group_name}/column_names", data=column_names, compression=compression)
                if index:
                    f.create_dataset(f"{group_name}/index", data=index, compression=compression)
        else:
            assert False, "Only pd.DataFrame and scipy.sparse.spmatrix are supported." + type(data)
        logger.info(f"Data {'appended' if mode=='a' else 'written'} to {path}")
    except Exception as e:
        logger.exception(e)
