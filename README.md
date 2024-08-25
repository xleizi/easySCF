
# README
中文版[这里](https://github.com/xleizi/easySCF/tree/main/README_zh.md)

This document describes how to perform file read and write operations in R and Python. In both R and Python, read and write operations are straightforward and can be accomplished with just one line of code.

## Read and Write in R

### Reading in R

In R, you can use the following function for file read and write operations:

- `loadH5()`: Used to read H5 files.

Parameter explanation:
- FileName: String, required. Specifies the path to the HDF5 file to be loaded.
- assay: String, default is "RNA". Specifies the type of layers to access in the HDF5 file.
- SeuratVersion: Provided by default by checkSeuratVersion(). This determines the version of the Seurat object used when loading data.
- image_name: String, default is "Spatial". Used to identify the name of the image data (if any).
- useBPcells: Boolean, default is FALSE. Whether to use BPcells for data storage.
- useDisk: Boolean, default is TRUE. Whether to use disk storage to optimize memory usage for BPcells data.
- calData: Boolean, default is TRUE. Whether to calculate normalization after loading data.
- calScale: Boolean, default is FALSE. Whether to scale the data.
- calFeatures: Boolean, default is FALSE. Whether to calculate highly variable genes.
- group_by: Default is NULL. Used to specify grouping variables, commonly used in subsequent analysis to differentiate between different data subsets.
- readType: String, default is "Seurat". Specifies the method or format of data reading, usually corresponding to the analysis software.

Here is a simple example demonstrating how to use these functions:

```R
library(easySCFr)
sce <- loadH5("data.h5")
```

### Writing in R

In R, you can use the following function for file read and write operations:

- `saveH5()`: Used to save Seurat objects to an H5 file.

Parameter explanation:
- data: Required. The Seurat object to be saved.
- FileName: String, required. Specifies the path and name of the file to save.
- assay: String, default is "RNA". Specifies the type of layers data contained in the file.
- save_graph: Boolean, default is TRUE. Whether to save data on cell interconnections.
- SeuratVersion: Provided by default by checkSeuratVersion(). This determines the version of the Seurat object used when saving data.
- image_name: String, default is NULL. If provided, specifies the name of the image file associated with the data.
- split_save: Boolean, default is TRUE. Whether to split the data into multiple subsets for saving, suitable for very large data sets.
- max_cells_per_subset: Integer, default is 5000. Specifies the maximum number of cells per subset during split saving.

Here is a simple example demonstrating how to use these functions:

```R
library(easySCFr)
saveH5(sce, "data.h5")
```

## Read and Write in Python

### Reading in Python

In Python, you can use the following function for file read and write operations:

- `loadH5()`: Used to read H5 files.

Parameter explanation:
- filename: str | Path, required. Specifies the path to the HDF5 file to be loaded.
- assay: str, default is "RNA". Specifies the type of layers data to be read from the HDF5 file.
- datatype: str, default is "scanpy". Specifies the data format, can be saved in different data types.
- image_name: str | None, optional. If provided, specifies the name of the image file associated with the data.
- backed: bool | Literal['r', 'r+'] | None, optional. Specifies the file read mode, 'r' for read-only mode, 'r+' for read-write mode, None for loading all into memory.
- as_sparse: Sequence[str], default is "raw/X". Specifies which data should be stored in a sparse matrix format.
- as_sparse_fmt: type[spmatrix], default is sparse.csr_matrix. Defines the format of the saved sparse matrix, modification is not recommended.

Here is a simple example demonstrating how to use these functions:

```python
from easySCF import loadH5
sce = loadH5("data.h5")
```

### Writing in Python

In Python, you can use the following function for file read and write operations:

- `saveH5()`: Used to save Seurat objects to an H5 file.

Parameter explanation:
- adata: Any, required. The AnnData object to be saved.
- h5_path: Path | str, required. The path for saving the HDF5 file.
- assay: str, default is "RNA". Specifies the type of layers data.
- datatype: str, default is "scanpy". Specifies the data format.
- image_name: str, default is "slice". Specifies the name of the image file associated with the data.
- save_graph: bool, default is True. Whether to save data on cell interconnections.
- as_dense: Sequence[str], default is an empty tuple. Manually specifies which data fields should be saved in a dense matrix format.
- split_save: bool, default is True. Whether to split the data into multiple subsets for saving, usually used for very large data sets.
- max_cells_per_subset: int, default is 5000. Specifies the maximum number of cells per subset.
- compression: Literal['gzip', 'lzf'] | None, default is "gzip". Specifies the compression algorithm used.
- compression_opts: int | None, optional. Provides specific parameters for the compression algorithm.

Here is a simple example demonstrating how to use these functions:

```python
from easySCF import saveH5
saveH5(sce, "data.h5")
```
