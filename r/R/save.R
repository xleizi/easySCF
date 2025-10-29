df_to_h5 <- function(h5data, df) {
  # 如果行数为 0,直接返回
  if (nrow(df) == 0) {
    return(NULL)
  }

  # 即使列数为 0,也保存 index (rownames)
  # 这对于 var/rawvar 非常重要,确保基因名称不会丢失
  h5AddAttribute(h5data, "encoding-type", "dataframe")
  h5AddAttribute(h5data, "shape", dim(df))

  # 如果有列,保存 column-order;否则保存空数组
  if (ncol(df) > 0) {
    h5AddAttribute(h5data, "column-order", colnames(df))
  } else {
    # 对于空 DataFrame,保存空数组作为 column-order
    h5AddAttribute(h5data, "column-order", character(0))
  }

  h5AddAttribute(h5data, "_index", "_index")
  h5AddAttribute(h5data, "encoding-version", "0.2.0")

  # 始终保存 index,即使 DataFrame 没有列
  h5dfAddIndex(h5data, df)

  # 如果没有列,直接返回
  if (ncol(df) == 0) {
    return(invisible(NULL))
  }

  # 保存每一列的数据
  for (name in colnames(df)) {
    if (name == "_index") {
      next
    }
    if (is.character(df[[name]]) & length(unique(df[[name]])) < 50) {
      df[[name]] <- as.factor(df[[name]])
    }
    if (is.factor(df[[name]])) {
      h5data2 <- h5data$create_group(name)
      df[[name]] <- droplevels(df[[name]])
      h5CreateStringDataset(h5data2, "categories", levels(df[[name]]))
      h5AddAttribute(h5data2[["categories"]], "encoding-type", "string-array")
      h5AddAttribute(h5data2[["categories"]], "encoding-version", "0.2.0")

      h5data2[["codes"]] <- as.integer(df[[name]]) - 1L
      h5AddAttribute(h5data2[["codes"]], "encoding-type", "array")
      h5AddAttribute(h5data2[["codes"]], "encoding-version", "0.2.0")

      h5AddAttribute(h5data2, "encoding-type", "categorical")
      h5AddAttribute(h5data2, "encoding-version", "0.2.0")
      h5AddAttribute(h5data2, "ordered", is.ordered(df[[name]]))
    } else {
      h5data[[name]] <- df[[name]]
      h5AddAttribute(h5data[[name]], "encoding-type", "array")
      h5AddAttribute(h5data[[name]], "encoding-version", "0.2.0")
    }
  }
}


X_to_h5 <- function(
    sce_X, h5data, data_name = "data",
    split_save = TRUE,
    max_cells_per_subset = 5000) {
  if (nrow(sce_X) == 0 | ncol(sce_X) == 0) {
    return(NULL)
  }
  if ("dgCMatrix" %in% class(sce_X)) {
    sce_X <- Matrix::t(sce_X)
    sce_X <- as(sce_X, "RsparseMatrix")
  }
  if ("Graph" %in% class(sce_X)) {
    sce_X <- as(sce_X, "RsparseMatrix")
  }
  # if ("matrix" %in% class(sce_X)) {
  #   # sce_X <- Matrix::t(sce_X)
  #   # sce_X <- sce_X
  # }
  if (!is.null(attr(class(sce_X), "package"))) {
    if (attr(class(sce_X), "package") == "BPCells") {
      sce_X <- Matrix::t(sce_X)
    }
  }

  handle_data_splitting(sce_X, h5data, data_name, split_save, max_cells_per_subset)
}

reductions_to_h5 <- function(sce, h5file, group_path) {
  reductionsh5 <- h5file$create_group(group_path)
  h5AddAttribute(reductionsh5, "encoding-type", "dict")
  h5AddAttribute(reductionsh5, "encoding-version", "0.1.0")
  for (name in names(sce@reductions)) {
    reductionsh5[[name]] <- t(sce@reductions[[name]]@cell.embeddings)
    h5AddAttribute(reductionsh5[[name]], "encoding-type", "array")
    h5AddAttribute(reductionsh5[[name]], "encoding-version", "0.2.0")
  }
}

graphs_to_h5 <- function(sce, h5file, group_path) {
  graphsh5 <- h5file$create_group(group_path)
  for (name in names(sce@graphs)) {
    tryCatch(
      {
        graphs_X <- as(sce@graphs[[name]], "RsparseMatrix")
        write_matrix(graphsh5, matrix = graphs_X, name = name)
      },
      error = function(e) {
        print(paste0("Error in writing graph: ", name, e))
      }
    )
  }
}

handle_data_splitting <- function(
    sce_X, h5data, data_name,
    split_save = TRUE,
    max_cells_per_subset = 5000) {
  if ("dgRMatrix" %in% class(sce_X)) {
    Xclass <- "csr_matrix"
  } else if ("dgCMatrix" %in% class(sce_X)) {
    Xclass <- "csc_matrix"
  } else if ("array" %in% class(sce_X)) {
    Xclass <- "array"
  } else {
    Xclass <- "csr_matrix"
  }

  h5data_Name <- h5data$create_group(data_name)
  h5AddAttribute(h5data_Name, "encoding-type", Xclass)
  h5AddAttribute(h5data_Name, "encoding-version", "0.1.0")
  h5AddAttribute(h5data_Name, "shape", dim(sce_X))

  num_cells <- nrow(sce_X)
  if (num_cells > 200000 & !split_save) {
    print("Large dataset detected. Splitting the save process.")
    split_save <- TRUE
  }
  if (split_save) {
    num_subsets <- ceiling(num_cells / max_cells_per_subset)
    # i <- 1
    for (i in 1:num_subsets) {
      start_idx <- ((i - 1) * max_cells_per_subset) + 1
      end_idx <- min(i * max_cells_per_subset, num_cells)
      subset <- sce_X[start_idx:end_idx, ]
      save_path <- sprintf("%s/%s_%05d.npz", data_name, data_name, i - 1)
      write_matrix(h5data, matrix = subset, name = save_path)
    }
  } else {
    write_matrix(h5data, matrix = sce_X, name = sprintf("%s/%s_00000.npz", data_name, data_name))
  }
}

# 插入稀疏矩阵
write_matrix <- function(h5data, matrix, name) {
  if (!is.null(attr(class(matrix), "package"))) {
    if (attr(class(matrix), "package") == "BPCells") {
      # matrix <- as(matrix, "dgRMatrix")
      # matrix <- as(matrix, "RsparseMatrix")
      matrix <- as(matrix, "sparseMatrix")
      matrix <- as(matrix, "RsparseMatrix")
    }
  }


  if ("matrix" %in% class(matrix)) {
    h5data[[name]] <- matrix
    h5AddAttribute(h5data[[name]], "encoding-type", "array")
    h5AddAttribute(h5data[[name]], "encoding-version", "0.2.0")
    h5AddAttribute(h5data[[name]], "shape", dim(matrix))
  } else if ("dgRMatrix" %in% class(matrix)) {
    h5data_sub <- h5data$create_group(name)
    h5data_sub[["data"]] <- slot(object = matrix, name = "x")
    h5data_sub[["indices"]] <- slot(object = matrix, name = "j")
    h5data_sub[["indptr"]] <- slot(object = matrix, name = "p")
    h5AddAttribute(h5data_sub, "encoding-type", "csr_matrix")
    h5AddAttribute(h5data_sub, "encoding-version", "0.1.0")
    h5AddAttribute(h5data_sub, "shape", dim(matrix))
  } else {
    print(paste0("Error in writing matrix: ", class(matrix)))
  }
}


# 将数据写入 HDF5 文件
commands_to_h5 <- function(data, h5file, group_path) {
  for (i in 1:length(data)) {
    if ("SeuratCommand" %in% class(data[[i]])) {
      tryCatch(
        {
          data[[i]] <- as.list(data[[i]])
        },
        error = function(e) {
          print(paste0("Error in converting SeuratCommand to list: ", e))
        }
      )
    }
  }
  # 创建或获取组
  if (!h5file$exists(group_path)) {
    group <- h5file$create_group(group_path)
  } else {
    group <- h5file$open_group(group_path)
  }

  # key <- names(data)[1]
  # 遍历列表数据
  for (key in names(data)) {
    item <- data[[key]]
    # 如果元素是列表，递归调用函数
    if (is.list(item)) {
      new_group_path <- paste0(group_path, "/", key)
      commands_to_h5(item, h5file, new_group_path)
    } else if (is.function(item)) {
      item <- capture.output(item)
      item <- paste(item, collapse = "\n")
      group[[key]] <- item
    } else if (is.function(item)) {
      data_type <- h5file$create_type("string", size = nchar(item))
      dataset <- group$create_dataset(key, dtype = data_type, dims = 1)
      group[[key]] <- item
    } else {
      group[[key]] <- item
    }
  }
}

#' saveH5
#'
#' @param FileName
#' @param data
#' @param assay
#' @param save_graph
#' @param SeuratVersion
#' @param split_save
#' @param max_cells_per_subset
#'
#' @return
#' @export
#'
#' @examples
saveH5 <- function(data, FileName, assay = "RNA",
                   save_graph = TRUE,
                   SeuratVersion = checkSeuratVersion(),
                   image_name = NULL,
                   split_save = TRUE,
                   max_cells_per_subset = 5000) {
  h5 <- hdf5r::H5File$new(FileName, "w")
  tryCatch(
    {
      data <- Seurat_to_H5(h5, data,
        assay = assay,
        save_graph = save_graph,
        SeuratVersion = SeuratVersion,
        split_save = split_save,
        image_name = image_name,
        max_cells_per_subset = max_cells_per_subset
      )
    },
    error = function(e) {
      print(e)
    },
    finally = {
      h5$close_all()
    }
  )

  if ("Seurat" %in% class(data)) {
    return(data)
  }
}

Seurat_to_H5 <- function(h5, sce, assay = "RNA",
                         save_graph = TRUE,
                         SeuratVersion = checkSeuratVersion(),
                         split_save = TRUE,
                         image_name = NULL,
                         max_cells_per_subset = 5000) {
  if (SeuratVersion == 5) {
    if ("Assay" %in% class(sce[[assay]])) {
      sce[[assay]] <- as(object = sce[[assay]], Class = "Assay5")
    }
    sce <- SeuratObject::JoinLayers(sce, assay = assay)
  }

  print("Saving X ...")
  # X
  assayList <- h5$create_group("assay")
  assayListsCur <- assayList$create_group(assay)
  layersList <- assayListsCur$create_group("layers")

  if (SeuratVersion == 4) {
    rawData <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "counts")
    scaleData <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "scale.data")
    # var
    varObjName <- "meta.features"
  } else if (SeuratVersion == 5) {
    rawData <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "counts")
    scaleData <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "scale.data")
    # var
    varObjName <- "meta.data"
  }
  X_to_h5(
    sce_X = rawData, h5data = layersList, data_name = "rawdata",
    split_save = split_save, max_cells_per_subset = max_cells_per_subset
  )

  X_to_h5(
    sce_X = scaleData, h5data = layersList, data_name = "data",
    split_save = split_save, max_cells_per_subset = max_cells_per_subset
  )

  print("Saving var ...")
  # var
  varList <- h5$create_group("var")

  rawvarh5 <- varList$create_group("rawvar")
  rawvar <- slot(object = sce@assays[[assay]], name = varObjName)

  # 如果 rawvar 为 NULL 或是空 DataFrame,创建一个只有 rownames 的空 DataFrame
  # 这确保即使没有 feature metadata,也能保存基因名称
  if (is.null(rawvar) || nrow(rawvar) == 0 || ncol(rawvar) == 0) {
    rawvar <- data.frame(row.names = rownames(sce@assays[[assay]]))
  } else {
    # 确保 rownames 正确
    rownames(rawvar) <- rownames(sce@assays[[assay]])
  }

  df_to_h5(rawvarh5, rawvar)

  # 检查是否有标准化数据（data 层）但没有 VariableFeatures
  hasScaleData <- (nrow(scaleData) > 0 && ncol(scaleData) > 0)

  if (!is.null(Seurat::VariableFeatures(sce))) {
    var <- rawvar[rownames(rawvar) %in% Seurat::VariableFeatures(sce), ]
    varh5 <- varList$create_group("var")
    df_to_h5(varh5, var)
  } else if (hasScaleData) {
    # 有标准化数据但没有 VariableFeatures，输出警告
    warning(
      "H5 file contains normalized/scaled data (data layer) but no VariableFeatures. ",
      "The 'var/var' metadata will not be saved, which may cause compatibility issues. ",
      "Consider running FindVariableFeatures() before saving for better compatibility.\n",
      "  Example: sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)"
    )
  }

  # obs
  print("Saving obs ...")
  meta <- sce@meta.data
  meta <- meta[colnames(rawData), ]
  for(col in colnames(meta)){
    if(any(is.na(meta[[col]]))){
        meta[[col]] <- "N/A"
    }
  }
  obs <- h5$create_group("obs")
  df_to_h5(obs, meta)

  # reductions
  print("Saving reductions ...")
  if (length(sce@reductions) > 0) {
    reductions_to_h5(sce, h5, "reductions")
    # tryCatch(
    #     {
    #         reductions_to_h5(sce, h5, "reductions")
    #     },
    #     error = function(e) {
    #         print(paste0("Error in writing reductions: ", e))
    #     }
    # )
  }

  # graphs
  print("Saving graphs ...")
  if (length(sce@graphs) > 0 & save_graph) {
    graphs_to_h5(sce, h5, "graphs")
  }

  # names
  h5CreateStringDataset(h5, "names_var", rownames(sce@assays[[assay]]))
  h5CreateStringDataset(h5, "names_obs", colnames(sce@assays[[assay]]))

  print("Saving commands ...")
  if (length(sce@commands) > 0) {
    commands_to_h5(sce@commands, h5, "uns")
  }

  print("Saving image ...")
  images_to_h5(sce, h5, image_name)
}


images_to_h5 <- function(sce, h5obj, image_name = NULL) {
  sceimages <- sce@images
  if(is.null(sceimages)){
    return(NULL)
  }
  if (length(sceimages) > 0) {
    if (!is.null(image_name)) {
      sceimage <- sceimages[[image_name]]
    } else {
      sceimage <- sceimages[[1]]
    }
  } else {
    return(NULL)
  }
  images <- h5obj$create_group("images")
  scale_factors <- sceimage@scale.factors
  images[["image"]] <- sceimage@image
  images[["coords"]] <- sceimage@boundaries[["centroids"]]@coords
  # images[["radius"]] <- sceimage@boundaries[["centroids"]]@radius
  scalefactors <- images$create_group("scale_factors")
  scalefactors[["spot"]] <- scale_factors$spot
  scalefactors[["fiducial"]] <- scale_factors$fiducial
  scalefactors[["hires"]] <- scale_factors$hires
  scalefactors[["lowres"]] <- scale_factors$lowres
}