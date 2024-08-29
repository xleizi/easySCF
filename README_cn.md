
# README

本文档介绍了如何在 R 和 Python 中进行文件读写操作。在 R 中和 Python 中的读写操作都很简单，均只需使用一句代码即可完成。

## 安装 easySCF

在 R 和 Python 中，你可以使用以下命令安装 easySCF：

```R
install.packages("devtools")
devtools::install_github("xleizi/easySCF/r")
```

```python
!pip install easySCFpy
```

## R 中的读写

### R 中的读取

在 R 中，你可以使用以下函数来进行文件的读写操作：

- `loadH5()`：用于读取 H5 文件。

参数说明:
- FileName: 字符串，必需。指定要加载的 HDF5 文件的路径。
- assay: 字符串，默认为 "RNA"。指定在 HDF5 文件中要访问的 layers 类型。
- SeuratVersion: 默认由 checkSeuratVersion() 提供。这决定了加载数据时使用的 Seurat 对象的版本。
- image_name: 字符串，默认为 "Spatial"。用于标识图像数据的名称（如果有的话）。
- useBPcells: 布尔值，默认为 FALSE。是否使用 BPcells 来储存数据。
- useDisk: 布尔值，默认为 TRUE。BPcells 数据是否使用磁盘存储来优化内存使用。
- calData: 布尔值，默认为 TRUE。是否在加载数据后计算 normalization。
- calScale: 布尔值，默认为 FALSE。是否对数据进行 scale。
- calFeatures: 布尔值，默认为 FALSE。是否计算高变基因。
- group_by: 默认为 NULL。用于指定分组变量，常用于后续分析中区分不同的数据子集。
- readType: 字符串，默认为 "Seurat"。指定读取数据的方式或格式，通常与分析软件对应。

以下是一个简单的示例，展示了如何使用这些函数：

```R
library(easySCFr)
sce <- loadH5("data.h5")
```

### R 中的写入

在 R 中，你可以使用以下函数来进行文件的读写操作：

- `saveH5()`：用于保存 Seurat 对象到 H5 文件。

参数说明:
- data: 必需。要保存的 Seurat 对象。
- FileName: 字符串，必需。指定保存文件的路径和名称。
- assay: 字符串，默认为 "RNA"。指定保存文件中包含的 layers 数据类型。
- save_graph: 布尔值，默认为 TRUE。是否保存细胞间联系的数据。
- SeuratVersion: 默认由 checkSeuratVersion() 提供。这决定了保存数据时使用的 Seurat 对象的版本。
- image_name: 字符串，默认为 NULL。如果提供，指定与数据关联的图像文件的名称。
- split_save: 布尔值，默认为 TRUE。是否将数据分割成多个子集来保存，适用于非常大的数据集。
- max_cells_per_subset: 整数，默认为 5000。在分割保存时，每个子集包含的最大细胞数。

以下是一个简单的示例，展示了如何使用这些函数：

```R
library(easySCFr)
saveH5(sce, "data.h5")
```

## Python 中的读写

### Python 中的读取

在 Python 中，你可以使用以下函数来进行文件的读写操作：

- `loadH5()`：用于读取 H5 文件。

参数说明:
- filename: str | Path，必需。指定要加载的 HDF5 文件的路径。
- assay: str，默认为 "RNA"。指定要从 HDF5 文件中读取的 layers 数据类型。
- datatype: str，默认为 "scanpy"。指定数据的格式，可保存为不同类型的数据。
- image_name: str | None，可选。如果提供，指定与数据关联的图像文件的名称。
- backed: bool | Literal['r', 'r+'] | None，可选。指定文件的读取模式，'r' 为只读模式，'r+' 为读写模式，None 表示全部载入内存。
- as_sparse: Sequence[str]，默认为 "raw/X"。指定哪些数据应以稀疏矩阵格式存储。
- as_sparse_fmt: type[spmatrix]，默认为 sparse.csr_matrix。定义保存的稀疏矩阵的格式，不建议修改。

以下是一个简单的示例，展示了如何使用这些函数：

```python
from easySCF import loadH5
sce = loadH5("data.h5")
```

### Python 中的写入

在 Python 中，你可以使用以下函数来进行文件的读写操作：

- `saveH5()`：用于保存 Seurat 对象到 H5 文件。

参数说明:
- adata: Any，必需。要保存的 AnnData 对象。
- h5_path: Path | str，必需。HDF5文件的保存路径。
- assay: str，默认为 "RNA"。指定layers数据类型。
- datatype: str，默认为 "scanpy"。指定数据的格式。
- image_name: str，默认为 "slice"。指定与数据关联的图像文件的名称。
- save_graph: bool，默认为 True。是否保存细胞间联系的数据。
- as_dense: Sequence[str]，默认为空元组。手动指定哪些数据字段应以密集矩阵格式保存。
- split_save: bool，默认为 True。是否将数据分割成多个子集进行保存，这通常用于非常大的数据集。
- max_cells_per_subset: int，默认为 5000。指定每个子集最多包含的细胞数。
- compression: Literal['gzip', 'lzf'] | None，默认为 "gzip"。指定使用的压缩算法。
- compression_opts: int | None，可选。提供压缩算法的具体参数。

以下是一个简单的示例，展示了如何使用这些函数：

```python
from easySCF import saveH5
saveH5(sce, "data.h5")
```
