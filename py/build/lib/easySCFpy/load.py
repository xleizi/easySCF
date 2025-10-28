import json
import numpy as np
from scipy import sparse
from anndata import AnnData
import h5py
from typing import Union, Sequence
from anndata._io.specs import read_elem
from anndata._io.h5ad import read_dataframe
from collections.abc import Sequence
from pathlib import Path
from typing import Literal


def read_dataframe_compat(h5_group):
    """
    兼容读取 DataFrame，支持旧版本文件中字符串格式的 column-order。

    :param h5_group: h5py.Group 对象，包含 DataFrame 数据
    :return: pandas DataFrame
    """
    import pandas as pd
    from anndata._io.specs import read_elem

    if not isinstance(h5_group, h5py.Group):
        raise TypeError("Expected h5py.Group")

    # 检查 column-order 属性
    col_order = h5_group.attrs.get('column-order')

    # 如果 column-order 是字符串（旧格式），手动读取
    if isinstance(col_order, str):
        columns = [col_order]
    elif isinstance(col_order, bytes):
        columns = [col_order.decode('utf-8')]
    elif isinstance(col_order, np.ndarray):
        # 正常格式，使用标准读取函数
        return read_dataframe(h5_group)
    else:
        # 没有 column-order 或其他格式，尝试标准读取
        return read_dataframe(h5_group)

    # 手动读取 DataFrame（处理字符串格式的 column-order）
    # 读取索引
    if '_index' in h5_group:
        index = read_elem(h5_group['_index'])
        if isinstance(index[0], bytes):
            index = [x.decode('utf-8') for x in index]
    else:
        index = None

    # 读取每一列数据
    data = {}
    for col in columns:
        if col in h5_group:
            data[col] = read_elem(h5_group[col])

    return pd.DataFrame(data, index=index)


def read_dense_matrix(h5_group):
    print("Reading dense matrix")
    denses = []
    for key in h5_group.keys():
        sub_matrix = read_elem(h5_group[key])
        denses.append(sub_matrix)

    # 如果只有一个矩阵，直接返回，避免不必要的 vstack
    if len(denses) == 1:
        return denses[0]

    if h5_group.attrs.get("encoding-type") == "csr_matrix":
        matrix = sparse.vstack(denses, format="csr")
    elif h5_group.attrs.get("encoding-type") == "csc_matrix":
        matrix = sparse.vstack(denses, format="csc")
    else:
        matrix = np.vstack(denses)
    return matrix


def read_matrice_matrix(h5_group, sparse_format):
    matrices = []
    for key in h5_group.keys():
        sub_matrix = read_elem(h5_group[key])
        matrices.append(sub_matrix)

    # 如果只有一个矩阵，直接返回，避免不必要的 vstack
    if len(matrices) == 1:
        return matrices[0]

    if sparse_format == sparse.csr_matrix:
        matrices = sparse.vstack(matrices, format="csr")
    else:
        matrices = sparse.vstack(matrices, format="csc")
    return matrices


def h5_to_X(layers, dataName, as_sparse, as_sparse_fmt, chunk_size):
    checkSp = "X" if dataName == "data" else "raw/X"
    if (
        not checkSp in as_sparse
        and layers[dataName].attrs.get("encoding-type") == "array"
    ):
        Data = read_dense_matrix(layers[dataName])
    else:
        Data = read_matrice_matrix(layers[dataName], sparse_format=as_sparse_fmt)
    return Data


def h5_to_misc(h5):
    data_dict = {}

    misc = h5["commands"]

    data_dict = {}

    def parse_group(group, parent_key=""):
        for key in group.keys():
            if isinstance(group[key], h5py.Dataset):
                data = group[key][()]
                data_dict[parent_key + key] = (
                    data if not isinstance(data, str) else data
                )
            elif isinstance(group[key], h5py.Group):
                parse_group(group[key], parent_key + key + "/")

    parse_group(misc)

    reorganized_data_dict = {}
    for key, value in data_dict.items():
        keys = key.split("/")
        current_dict = reorganized_data_dict
        for k in keys[:-1]:
            current_dict = current_dict.setdefault(k, {})
        current_dict[keys[-1]] = value

    def deserialize_dict(d):
        for k, v in d.items():
            if isinstance(v, dict):
                d[k] = deserialize_dict(v)
            elif isinstance(v, str):
                try:
                    d[k] = json.loads(v)
                except ValueError:
                    pass
        return d

    reorganized_data_dict = deserialize_dict(reorganized_data_dict)

    return reorganized_data_dict


def h5_to_uns_dict(h5file: h5py.File, path: str) -> dict:
    """
    Recursively reads a dictionary from an h5 file.

    :param h5file: h5py file object
    :param path: path in the h5 file where the data is stored
    :return: dictionary read from the h5 file
    """
    data = {}
    for key in h5file[path].keys():
        new_path = f"{path}/{key}"
        if isinstance(h5file[new_path], h5py.Group):
            data[key] = h5_to_uns_dict(h5file, new_path)
        else:
            item = h5file[new_path][()]
            if isinstance(item, bytes):
                try:
                    item = json.loads(item.decode("utf-8"))
                except json.JSONDecodeError:
                    item = item.decode("utf-8")
            data[key] = item
    return data


def read_h5_to_scanpy(
    filename: Union[str, Path],
    assay: str = "RNA",
    image_name: str | None = None,
    backed: Union[Literal["r"], Literal["r+"], bool, None] = None,
    *,
    as_sparse: Sequence[str] = ("raw/X"),
    as_sparse_fmt: type[sparse.spmatrix] = sparse.csr_matrix,
    chunk_size: int = 6000,
) -> AnnData:

    if as_sparse_fmt not in (sparse.csr_matrix, sparse.csc_matrix):
        raise NotImplementedError(
            "Dense formats can only be read to CSR or CSC matrices at this time."
        )
    if isinstance(as_sparse, str):
        as_sparse = [as_sparse]
    else:
        as_sparse = list(as_sparse)
    for i in range(len(as_sparse)):
        if as_sparse[i] in {("raw", "X"), "raw.X"}:
            as_sparse[i] = "raw/X"
        elif as_sparse[i] not in {"raw/X", "X"}:
            raise NotImplementedError(
                "Currently only `X` and `raw/X` can be read as sparse."
            )
    with h5py.File(filename, "r") as h5:
        if "var" not in h5.keys():
            raise KeyError("var not found.")
        if "obs" not in h5.keys():
            raise KeyError("obs not found.")
        if "assay" not in h5.keys():
            raise KeyError("assay not found.")
        try:
            layers = h5.get("assay", {}).get(assay, {}).get("layers", None)
            print("Reading X data")
            if layers is not None:
                if "data" in layers.keys():
                    data = h5_to_X(layers, "data", as_sparse, as_sparse_fmt, chunk_size)
                    rawData = h5_to_X(
                        layers, "rawdata", as_sparse, as_sparse_fmt, chunk_size
                    )
                    obs = read_dataframe_compat(h5["obs"])
                    obs = obs if obs.shape[0] != 0 and obs.shape[1] != 0 else None
                    # 检查 var/var 是否存在，不存在则回退到 var/rawvar
                    if "var" in h5["var"].keys():
                        var = read_dataframe_compat(h5["var/var"])
                    else:
                        var = read_dataframe_compat(h5["var/rawvar"])
                    var = var if var.shape[0] != 0 and var.shape[1] != 0 else None
                    adata = AnnData(
                        X=data,
                        obs=obs,
                        var=var,
                    )

                    obs = read_dataframe_compat(h5["obs"])
                    obs = obs if obs.shape[0] != 0 and obs.shape[1] != 0 else None
                    var = read_dataframe_compat(h5["var/rawvar"])
                    var = var if var.shape[0] != 0 and var.shape[1] != 0 else None
                    adata_raw = AnnData(
                        X=rawData,
                        obs=obs,
                        var=var,
                    )
                    adata.raw = adata_raw
                else:
                    rawData = h5_to_X(
                        layers, "rawdata", as_sparse, as_sparse_fmt, chunk_size
                    )
                    adata = AnnData(
                        X=rawData,
                        obs=read_dataframe_compat(h5["obs"]),
                        var=read_dataframe_compat(h5["var/rawvar"]),
                    )
            else:
                raise KeyError("Layers not found.")
        except Exception as e:
            import traceback
            print(f"详细错误信息: {e}")
            print(f"错误类型: {type(e).__name__}")
            traceback.print_exc()
            raise e
        print("Reading reductions data")
        if "reductions" in h5.keys():
            adata.obsm = read_elem(h5["reductions"])
        print("Reading graphs data")
        if "graphs" in h5.keys():
            adata.obsp = read_elem(h5["graphs"])
        print("Reading commands data")
        if "commands" in h5.keys():
            adata.uns = h5_to_uns_dict(h5, "uns")
        print("Reading images data")
        if "images" in h5.keys():
            adata = h5_to_image(adata, h5["images"], image_name)
    return adata


def loadH5(
    filename: Union[str, Path],
    assay: str = "RNA",
    datatype: str = "scanpy",
    image_name: str | None = None,
    backed: Union[Literal["r"], Literal["r+"], bool, None] = None,
    *,
    as_sparse: Sequence[str] = ("raw/X"),
    as_sparse_fmt: type[sparse.spmatrix] = sparse.csr_matrix,
    chunk_size: int = 6000,
) -> AnnData:
    if datatype == "scanpy":
        return read_h5_to_scanpy(
            filename,
            assay=assay,
            image_name=image_name,
            backed=backed,
            as_sparse=as_sparse,
            as_sparse_fmt=as_sparse_fmt,
            chunk_size=chunk_size,
        )
    else:
        raise NotImplementedError(f"Datatype {datatype} not supported.")


def h5_to_image(adata: AnnData, image_layers: h5py.Group, image_name: str | None = None) -> AnnData:
    if image_name is None:
        image_name = "slice"
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][image_name] = dict()
    adata.uns["spatial"][image_name]["images"] = dict()

    image = np.transpose(read_elem(image_layers["image"]), (2, 1, 0))
    # image = np.transpose(image_layers["image"][:], (1, 2, 0))
    if len(image[1]) > 2000:
        res = "hires"
    else:
        res = "lowres"
    adata.uns["spatial"][image_name]["images"][res] = image

    fiducial = read_elem(image_layers["scale_factors/fiducial"])
    hires = read_elem(image_layers["scale_factors/hires"])
    lowres = read_elem(image_layers["scale_factors/lowres"])
    spot = read_elem(image_layers["scale_factors/spot"])
    
    if isinstance(fiducial, np.ndarray):
        fiducial = fiducial[0]
    if isinstance(hires, np.ndarray):
        hires = hires[0]
    if isinstance(lowres, np.ndarray):
        lowres = lowres[0]
    if isinstance(spot, np.ndarray):
        spot = spot[0]

    scalefactors = {
        "spot_diameter_fullres": spot,
        # "bin_size_um": 8.0,
        # "microns_per_pixel": 0.27389195078106876,
        # "regist_target_img_scalef": 0.07973422,
        "tissue_lowres_scalef": lowres,
        "fiducial_diameter_fullres": fiducial,
        "tissue_hires_scalef": hires,
    }

    adata.uns["spatial"][image_name]["scalefactors"] = scalefactors
    transposed_coords = np.transpose(read_elem(image_layers["coords"]), (1, 0))
    transposed_coords[:, [0, 1]] = transposed_coords[:, [1, 0]]
    adata.obsm["spatial"] = transposed_coords
    return adata
