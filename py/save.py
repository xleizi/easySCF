import json
import h5py
from anndata import AnnData
from typing import Literal, Optional, Union
from pathlib import Path
import numpy as np
from scipy import sparse
from anndata._io.specs import write_elem
from anndata._core.sparse_dataset import BaseCompressedSparseDataset
from collections.abc import Mapping, Sequence
from types import MappingProxyType
import scanpy as sc


def json_serialize(obj):
    if isinstance(obj, np.int64):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def uns_to_h5(
    h5file: h5py.File,
    data: dict,
    path: str,
):
    """
    Recursively writes a dictionary to an h5 file.

    :param h5file: h5py file object
    :param path: path in the h5 file where the data should be written
    :param data: dictionary to write to the h5 file
    """
    for key, item in data.items():
        new_path = f"{path}/{key}"
        if isinstance(item, dict):
            uns_to_h5(h5file, item, new_path)
        else:
            if isinstance(item, bytes):
                str_data = item.decode("utf-8")
                tmp = json.loads(str_data)
                tmp_json = json.dumps(tmp, default=json_serialize)
                h5file.create_dataset(new_path, data=tmp_json)
            elif isinstance(
                item,
                (np.ndarray, np.generic, int, float),
            ):
                h5file.create_dataset(new_path, data=item)
            elif isinstance(item, str):
                h5file.create_dataset(new_path, data=item.encode("utf-8"))
            else:
                h5file.create_dataset(new_path, data=str(item))


def obs_to_h5(
    adata: AnnData,
    h5: h5py.File,
    obsname="obs",
    dataset_kwargs: Mapping = MappingProxyType({}),
):
    if hasattr(adata, obsname):
        adata.obs_names.name = None
        write_elem(h5, obsname, adata.obs, dataset_kwargs=dataset_kwargs)


def obsm_to_h5(
    adata: AnnData,
    h5: h5py.File,
    obsmname="obsm",
    dataset_kwargs: Mapping = MappingProxyType({}),
):

    if hasattr(adata, obsmname):
        for i in adata.obsm:
            if type(adata.obsm[i]) == "pandas.core.frame.DataFrame":
                adata.obsm[i] = adata.obsm[i].values
        write_elem(h5, "reductions", dict(adata.obsm), dataset_kwargs=dataset_kwargs)


def obsp_to_h5(
    adata: AnnData,
    h5: h5py.File,
    save_graph,
    obspname="obsp",
    dataset_kwargs: Mapping = MappingProxyType({}),
):
    if save_graph and hasattr(adata, obspname):
        write_elem(h5, "graphs", dict(adata.obsp), dataset_kwargs=dataset_kwargs)


def spatial_to_h5(
    adata: AnnData,
    h5: h5py.File,
    spatialname="spatial",
    dataset_kwargs: Mapping = MappingProxyType({}),
):
    if hasattr(adata, spatialname):
        write_elem(h5, spatialname, dict(adata.spatial), dataset_kwargs=dataset_kwargs)


def handle_data_splitting(
    adata_X: AnnData,
    layers: h5py.Group,
    data_name: str,
    split_save: bool,
    max_cells_per_subset: int,
    dataset_kwargs: Mapping,
):
    num_cells = adata_X.shape[0]
    if num_cells > 200000 and not split_save:
        print("Large dataset detected. Splitting the save process.")
        split_save = True

    if split_save:
        num_subsets = np.ceil(num_cells / max_cells_per_subset).astype(int)
        for i in range(num_subsets):
            start_idx = i * max_cells_per_subset
            end_idx = min((i + 1) * max_cells_per_subset, num_cells)
            subset = adata_X[start_idx:end_idx]
            save_path = f"{data_name}/{data_name}_{i:05d}.npz"
            write_elem(layers, save_path, subset, dataset_kwargs=dataset_kwargs)
    else:
        write_elem(
            layers,
            f"{data_name}/{data_name}_00000.npz",
            adata_X,
            dataset_kwargs=dataset_kwargs,
        )
    layers[data_name].attrs.update(
        {
            "encoding-type": layers[f"{data_name}/{data_name}_00000.npz"].attrs[
                "encoding-type"
            ],
            "encoding-version": layers[f"{data_name}/{data_name}_00000.npz"].attrs[
                "encoding-version"
            ],
            "shape": adata_X.shape,
        }
    )


def write_X(
    adata_X: AnnData,
    layers: h5py.Group,
    X_name: str,
    data_name: str,
    as_dense: Sequence[str] = (),  # 是否存为密集矩阵
    split_save: bool = True,
    max_cells_per_subset=5000,
    dataset_kwargs: Mapping = MappingProxyType({}),
):
    if X_name in as_dense and isinstance(
        adata_X, (sparse.spmatrix, BaseCompressedSparseDataset)
    ):
        try:
            adata_X = adata_X.toarray()
        except AttributeError:
            pass

    handle_data_splitting(
        adata_X, layers, data_name, split_save, max_cells_per_subset, dataset_kwargs
    )


def X_to_h5(
    adata: AnnData,
    layers: h5py.Group,
    var: h5py.Group,
    as_dense: Sequence[str] = (),  # 是否存为密集矩阵
    splite_save: bool = True,
    max_cells_per_subset=5000,
    dataset_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
):
    dataset_kwargs = {**dataset_kwargs, **kwargs}
    if adata.raw is not None:
        write_X(
            adata.X,
            layers,
            "X",
            "data",
            as_dense,
            splite_save,
            max_cells_per_subset,
            dataset_kwargs,
        )
        write_X(
            adata.raw.X,
            layers,
            "raw/X",
            "rawdata",
            as_dense,
            splite_save,
            max_cells_per_subset,
            dataset_kwargs,
        )

        # 存储var和raw var
        adata.var_names.name = None
        adata.raw.var_names.name = None
        write_elem(var, "var", adata.var, dataset_kwargs=dataset_kwargs)
        write_elem(var, "rawvar", adata.raw.var, dataset_kwargs=dataset_kwargs)
    else:
        write_X(
            adata.X,
            layers,
            "X",
            "rawdata",
            as_dense,
            splite_save,
            max_cells_per_subset,
            dataset_kwargs,
        )

        # 存储var
        adata.var_names.name = None
        write_elem(var, "rawvar", adata.var, dataset_kwargs=dataset_kwargs)


def saveH5(
    adata,
    h5_path: Union[Path, str],
    assay: str = "RNA",
    datatype: str = "scanpy",
    image_name: str = "slice",
    save_graph: bool = True,
    as_dense: Sequence[str] = (),
    split_save: bool = True,
    max_cells_per_subset: int = 5000,
    compression: Optional[Literal["gzip", "lzf"]] = "gzip",
    compression_opts: Optional[int] = None,
    dataset_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
):
    if datatype == "scanpy":
        scanpy_write_to_h5(
            adata,
            h5_path,
            assay,
            image_name, 
            save_graph,
            as_dense,
            split_save,
            max_cells_per_subset,
            compression,
            compression_opts,
            dataset_kwargs,
            **kwargs,
        )


def scanpy_write_to_h5(
    adata: AnnData,
    h5_path: Union[Path, str],
    assay: str = "RNA",
    image_name: str = "slice",
    save_graph: bool = True,
    as_dense: Sequence[str] = (),
    split_save: bool = True,
    max_cells_per_subset: int = 5000,
    compression: Optional[Literal["gzip", "lzf"]] = "gzip",
    compression_opts: Optional[int] = None,
    dataset_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
):
    if isinstance(as_dense, str):
        as_dense = [as_dense]
    if "raw.X" in as_dense:
        as_dense = list(as_dense)
        as_dense[as_dense.index("raw.X")] = "raw/X"
    if any(val not in {"X", "raw/X"} for val in as_dense):
        raise NotImplementedError(
            "Currently, only `X` and `raw/X` are supported values in `as_dense`"
        )
    if "raw/X" in as_dense and adata.raw is None:
        raise ValueError("Cannot specify writing `raw/X` to dense if it doesn’t exist.")

    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = {**dataset_kwargs, **kwargs}
    h5_path = Path(h5_path)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()

    # 将路径字符串转换为Path对象，以确保兼容性
    with h5py.File(h5_path, mode) as h5:
        h5.attrs["app"] = "scanpy"
        h5.attrs["h5app"] = "h5py"
        h5.attrs[h5.attrs["app"]] = sc.__version__
        h5.attrs[h5.attrs["h5app"]] = h5py.__version__
        layers = h5.create_group(f"assay/{assay}/layers")

        var = h5.create_group("var")

        print("Saving X ...")
        X_to_h5(
            adata,
            layers,
            var,
            as_dense,
            split_save,
            max_cells_per_subset,
            dataset_kwargs,
            compression=compression,
            compression_opts=compression_opts,
        )

        print("Saving obs ...")
        # 保存obs
        obs_to_h5(adata, h5, "obs", dataset_kwargs)

        print("Saving spatial ...")
        if "spatial" in adata.uns and "spatial" in adata.obsm:
            images = h5.create_group("images")
            image_to_h5(images, adata, image_name)
            del adata.obsm["spatial"]
            del adata.uns["spatial"]

        # dimR即reduction部分obsm
        # 保存obsm
        print("Saving obsm ...")
        obsm_to_h5(
            adata,
            h5,
            "obsm",
            dataset_kwargs,
        )

        print("Saving obsp ...")
        obsp_to_h5(adata, h5, save_graph, "obsp", dataset_kwargs)

        h5.create_dataset("names_var", data=adata.var_names, dtype=h5py.string_dtype())
        h5.create_dataset("names_obs", data=adata.obs_names, dtype=h5py.string_dtype())
        print("Saving uns ...")
        # 保存uns
        uns_to_h5(
            h5,
            adata.uns,
            "uns",
        )


def image_to_h5(layers: h5py.Group, adata: AnnData, image_name: str = "slice"):
    imageList = adata.uns["spatial"][image_name]["images"]
    if len(imageList) > 1:
        print("More than one image found, but only one image will be used")
    if "lowres" in imageList.keys():
        image = imageList["lowres"]
        print("Using lowres image")
    elif "hires" in imageList.keys():
        image = imageList["hires"]
        print("Using hires image")
    else:
        print("No lowres or hires image found, using the first one")
        image = imageList[next(iter(imageList))]
    image = np.transpose(image, (2, 1, 0))

    coords = adata.obsm["spatial"]
    coords[:, [0, 1]] = coords[:, [1, 0]]
    coords = np.transpose(coords, (1, 0))

    write_elem(layers, "image", image)
    write_elem(layers, "coords", coords)

    scale_factorsGroup = layers.create_group("scale_factors")

    scale_factors = adata.uns["spatial"][image_name]["scalefactors"]
    spot = scale_factors["spot_diameter_fullres"]
    fiducial = scale_factors["fiducial_diameter_fullres"]
    hires = scale_factors["tissue_hires_scalef"]
    lowres = scale_factors["tissue_lowres_scalef"]
    write_elem(scale_factorsGroup, "spot", spot)
    write_elem(scale_factorsGroup, "fiducial", fiducial)
    write_elem(scale_factorsGroup, "hires", hires)
    write_elem(scale_factorsGroup, "lowres", lowres)
