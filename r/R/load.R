read_matrice_matrix <- function(h5data) {
  if ("data" %in% names(h5data) && "indices" %in% names(h5data) && "indptr" %in% names(h5data)) {
    i <- h5data[["indices"]][]
    p <- h5data[["indptr"]][]
    x <- h5data[["data"]][]

    X <- if ("csr_matrix" %in% getEncodingType(h5data)) {
      Matrix::sparseMatrix(i = i, p = p, x = x, dims = rev(getShape(h5data)), index1 = FALSE)
    } else {
      Matrix::t(Matrix::sparseMatrix(j = i, p = p, x = x, dims = getShape(h5data), index1 = FALSE))
    }

    return(X)
  } else {
    return(h5data[])
  }
}

# name <- names(h5layer)[[1]]
# h5data <- group_h5layer
read_matrix <- function(h5layer, useBPcells = FALSE, cellNames = NULL, geneNames = NULL) {
  data_list <- list()
  for (name in names(h5layer)) {
    group_h5layer <- h5layer[[name]]
    if ("array" %in% getEncodingType(h5layer)) {
      data_list[[name]] <- group_h5layer[, ]
    } else {
      if (useBPcells) {
        data_list[[name]] <- as(read_matrice_matrix(group_h5layer), "IterableMatrix")
      } else {
        data_list[[name]] <- read_matrice_matrix(h5data = group_h5layer)
      }
    }
  }

  # dim(data_list[[1]])
  # dim(data_list[[2]])
  # rbind(data_list[[1]], data_list[[2]])
  all_data <- do.call(cbind, data_list)
  if (!is.null(cellNames) & ncol(all_data) == length(cellNames)) {
    colnames(all_data) <- cellNames
  }
  if (!is.null(geneNames) & nrow(all_data) == length(geneNames)) {
    rownames(all_data) <- geneNames
  }
  return(all_data)
}


h5_to_X <- function(h5, assay = "RNA", layer = "rawdata", useBPcells = FALSE, useDisk = TRUE, cellNames = NULL, geneNames = NULL) {
  if (is.null(cellNames)) {
    cellNames <- h5[["names_obs"]][]
  }

  if (is.null(geneNames)) {
    if (layer == "rawdata") {
      geneNames <- h5[["var/rawvar/_index"]][]
    } else {
      geneNames <- h5[["names_var"]][]
    }
  }

  h5layer <- h5[["assay"]][[assay]][["layers"]][[layer]]
  if (length(names(h5layer)) > 200) {
    print("The number of layers is too large. using BPCells to load the data.")
    useBPcells <- TRUE
  }
  if (useBPcells) {
    if (!require(BPCells)) {
      stop("The number of layers is too large. But the BPCells package not installed.")
    }
  }

  varType <- if (layer == "rawdata") "rawvar" else "var"
  X <- read_matrix(h5layer, useBPcells = useBPcells, cellNames = cellNames, geneNames = geneNames)
  if (useDisk & useBPcells) {
    diskFile <- sprintf("./.BPCells_dont_delete/%s/%s/%s", sub(".h5", "", basename(h5$filename)), assay, layer)
    if (dir.exists(diskFile)) {
      print("The data has been cached in disk. Delete the cache file to resave the data.")
      unlink(diskFile, recursive = TRUE)
    }
    write_matrix_dir(X, diskFile)
    X <- open_matrix_dir(diskFile)
  }
  return(X)
}

h5_to_DF <- function(h5data) {
  if (h5data[["_index"]]$dims == 0) {
    return(NULL)
  }
  rownamesStr <- NULL
  for (name in names(h5data)) {
    if (name == "_index") {
      rownamesStr <- h5data[[name]][]
      next
    }
  }
  df <- data.frame()
  for (name in names(h5data)) {
    if (name == "_index") {
      next
    }
    if ("categorical" %in% getEncodingType(h5data[[name]])) {
      values_attr <- h5data[[name]]
      labelName <- values_attr[["categories"]][]
      values <- values_attr[["codes"]][]
      values <- factor(as.integer(values), labels = labelName)
    } else if (getEncodingType(h5data[[name]]) %in% c("array", "string-array")) {
      values <- h5data[[name]][]
    } else {
      stop("Unknown encoding type")
    }
    df <- addDF(df, setNames(data.frame(values), name), "col")
  }
  if (!is.null(rownamesStr)) {
    rownames(df) <- rownamesStr
  }
  if ("column-order" %in% hdf5r::h5attr_names(h5data)) {
    colnamesOrder <- hdf5r::h5attr(h5data, "column-order")
    df <- df[, colnamesOrder]
  }

  return(df)
}


openH5 <- function(FileName) {
  if (!hdf5r::is_hdf5(FileName)) {
    stop("File is not a hdf5 file.")
  }
  h5 <- hdf5r::H5File$new(FileName, "r")
  if (!("var" %in% names(h5))) {
    stop("var not found.")
  }
  if (!("obs" %in% names(h5))) {
    stop("obs not found.")
  }
  if (!("assay" %in% names(h5))) {
    stop("assay not found.")
  }
  return(h5)
}

# 读取 uns 数据
h5_to_uns <- function(h5obj) {
  data <- list()

  # 检查 h5obj 是否为组
  if (is(h5obj, "H5Group")) {
    keys <- names(h5obj)
    for (key in keys) {
      # 使用 tryCatch 跳过报错项
      tryCatch(
        {
          data[[key]] <- h5_to_uns(h5obj[[key]])
        },
        error = function(e) {
          message(paste("message: Error processing key:", key, "-", e$message))
        }
      )
    }
  } else { # h5obj 是数据集
    item <- tryCatch(
      {
        h5obj[]
      },
      error = function(e) {
        message(paste("message: Error reading dataset:", e$message))
        NULL
      }
    )

    if (!is.null(item)) {
      if (is.raw(item[1])) {
        item_text <- rawToChar(item)
        # 尝试将内容作为 JSON 解析，否则保持原样
        data <- tryCatch(
          {
            jsonlite::fromJSON(item_text)
          },
          error = function(e) {
            message("message: Error parsing JSON - returning raw text")
            item_text
          }
        )
      } else {
        data <- item
      }
    }
  }
  return(data)
}


h5_to_obs <- function(h5, obsName = "obs") {
  return(h5_to_DF(h5[["obs"]]))
}

sce_add_h5_to_graphs <- function(sce, h5, cellNames, graphsName = "graphs") {
  for (name in names(h5[[graphsName]])) {
    tryCatch(
      {
        Graphmt <- read_matrice_matrix(h5[[graphsName]][[name]])
        rownames(Graphmt) <- cellNames
        colnames(Graphmt) <- cellNames
        sce@graphs[[name]] <- Seurat::as.Graph(Graphmt)
      },
      error = function(e) {
        print(paste0("Reading graph ", name, ": ", e$message, "failed."))
      }
    )
  }
  return(sce)
}

sce_add_h5_to_var <- function(sce, h5, assay, varName = "var", SeuratVersion = checkSeuratVersion()) {
  if (SeuratVersion == 5) {
    sce@assays[[assay]] <- Seurat::AddMetaData(sce@assays[[assay]], h5_to_DF(h5[[varName]][["rawvar"]]))
  } else if (SeuratVersion == 4) {
    sce@assays[[assay]]@meta.features <- h5_to_DF(h5[[varName]][["rawvar"]])
  } else {
    stop("Unsupported Seurat version")
  }
  return(sce)
}

sce_add_h5_to_reductions <- function(sce, h5, cellNames, assay = "RNA", reductionsName = "reductions") {
  name <- names(h5[[reductionsName]])[[1]]
  for (name in names(h5[[reductionsName]])) {
    tryCatch(
      {
        matrix <- h5[[reductionsName]][[name]][, ]
        colnames(matrix) <- cellNames
        sce@reductions[[name]] <- Seurat::CreateDimReducObject(
          embeddings = t(matrix),
          key = paste0(name, "_"),
          assay = assay
        )
      },
      error = function(e) {
        print(paste0("Reading reductions ", name, ": ", e$message, "failed."))
      }
    )
  }
  return(sce)
}

#' loadH5
#'
#' @param FileName
#' @param assay
#' @param SeuratVersion
#' @param calData
#' @param calScale
#' @param calFeatures
#' @param readType
#'
#' @return
#' @export
#'
#' @examples
loadH5 <- function(FileName,
                   assay = "RNA",
                   SeuratVersion = checkSeuratVersion(),
                   image_name = "Spatial",
                   useBPcells = FALSE,
                   useDisk = TRUE,
                   calData = TRUE,
                   calScale = FALSE,
                   calFeatures = FALSE,
                   group_by = NULL,
                   readType = "Seurat") {
  options(warn = -1)
  if (!readType %in% c("Seurat", "monocle2", "monocle3", "cellchat", "SingleCellExperiment")) {
    stop("readType should be one of Seurat, monocle2, monocle3")
  }

  if (readType == "monocle2") {
    if (!require(monocle)) {
      stop("The monocle package is not installed, please install it first.")
    }
  }

  if (readType == "monocle3") {
    if (!require(monocle3)) {
      stop("The monocle3 package is not installed, please install it first.")
    }
  }

  if (readType == "cellchat") {
    if (!require(CellChat)) {
      stop("The cellchat package is not installed, please install it first.")
    }
    if (is.null(group_by)) {
      stop("group_by should be specified for cellchat data.")
    }
  }

  if (readType == "SingleCellExperiment") {
    if (!require(SingleCellExperiment)) {
      stop("The SingleCellExperiment package is not installed, please install it first.")
    }
  }

  h5 <- openH5(FileName)
  tryCatch(
    {
      sce <- h5_to_seurat(h5, assay, SeuratVersion,
        image_name = image_name,
        useBPcells = useBPcells, useDisk = useDisk, calData = calData, calScale = calScale, calFeatures = calFeatures
      )
    },
    error = function(e) {
      print(e)
    },
    finally = {
      h5$close_all()
    }
  )

  if (readType == "Seurat") {
    return(sce)
  } else if (readType == "monocle2") {
    library(monocle)
    cds <- Seurat_to_Monocle2(sce, assay = assay, SeuratVersion = SeuratVersion)
    return(cds)
  } else if (readType == "monocle3") {
    library(monocle3)
    cds <- Seurat_to_Monocle3(sce, assay = assay, SeuratVersion = SeuratVersion)
    return(cds)
  } else if (readType == "cellchat") {
    loadlayers <- "data"
    if (SeuratVersion == 5) {
      if (!("data" %in% names(sce@assays[[assay]]@layers)) & !calData) {
        print("No data found in the h5 file and calData is FALSE, using the counts in the SeuratObject.")
        loadlayers <- "counts"
      }
    }
    cellchatdata <- Seurat_to_cellchat(sce, assay = assay, SeuratVersion = SeuratVersion, group_by = group_by, loadlayers = loadlayers)
    return(cellchatdata)
  } else if (readType == "SingleCellExperiment") {
    SingleCellData <- Seurat::as.SingleCellExperiment(sce)
    return(SingleCellData)
  } else {
    return(sce)
  }
}

Seurat_to_cellchat <- function(
    sce,
    assay = "RNA",
    SeuratVersion = checkSeuratVersion(),
    group_by = "orig.ident",
    loadlayers = "data") {
  if (SeuratVersion == 4) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = loadlayers)
  } else if (SeuratVersion == 5) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = loadlayers)
  }
  meta <- sce@meta.data
  cellchat <- createCellChat(object = data, meta = meta, group.by = group_by)
  return(cellchat)
}

Seurat_to_Monocle3 <- function(
    sce,
    assay = "RNA",
    SeuratVersion = checkSeuratVersion()) {
  if (SeuratVersion == 4) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "counts")
  } else if (SeuratVersion == 5) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "counts")
  }

  cell_metadata <- sce@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(sce))
  rownames(gene_annotation) <- rownames(sce)

  cds <- monocle3::new_cell_data_set(data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  return(cds)
}

Seurat_to_Monocle2 <- function(
    sce,
    assay = "RNA",
    SeuratVersion = checkSeuratVersion()) {
  if (SeuratVersion == 4) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "counts")
  } else if (SeuratVersion == 5) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "counts")
  }

  data <- as(as.matrix(data), "sparseMatrix")
  pd <- new("AnnotatedDataFrame", data = sce@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)

  cds <- monocle::newCellDataSet(data,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.1,
    expressionFamily = negbinomial.size()
  )
  return(cds)
}

h5_to_seurat <- function(h5,
                         assay = "RNA",
                         SeuratVersion = checkSeuratVersion(),
                         image_name = "Spatial",
                         useBPcells = FALSE,
                         useDisk = TRUE,
                         calData = TRUE,
                         calScale = FALSE,
                         calFeatures = FALSE) {
  cellNames <- h5[["names_obs"]][]
  geneNames <- h5[["names_var"]][]
  allgeneNames <- h5[["var/rawvar/_index"]][]
  print(sprintf("there is assay: %s in the h5 file", names(h5[["assay"]])))
  # 1 如果存在data, data ->"" sce@data, rawdata -> sce@counts,var的行名 -> sce@assays[["RNA"]]@var.feature
  # 2 如果不存在data, rawdata -> sce@counts
  # assay
  layersNames <- names(h5[["assay"]][[assay]][["layers"]])
  print("Reading raw data...")
  # layer <- "rawdata"
  rawX <- h5_to_X(h5, assay, layer = "rawdata", useBPcells = useBPcells, useDisk = useDisk, cellNames = cellNames, geneNames = allgeneNames)
  sce <- Seurat::CreateSeuratObject(counts = rawX, assay = assay, project = "scanpy")
  print("Normalization counts...")
  if (calData) {
    sce <- Seurat::NormalizeData(sce)
  }

  print("Add var...")
  # var
  sce <- sce_add_h5_to_var(sce, h5, assay, "var")

  print("Add Feature...")
  # Features
  if (calFeatures) {
    # if (is.null(sce@assays[[assay]]@layers$scale.data)) {
    #   if ("data" %in% layersNames) {
    #     print("there is no scale data in h5 files and calScale is FALSE, skip FindVariableFeatures")
    #     Seurat::VariableFeatures(sce) <- rownames(h5_to_DF(h5[["var"]][["var"]]))
    #   } else {
    #     print("there is no scale data in h5 files and calScale is FALSE, please set both 'calData = TRUE' and 'calScale = TRUE' to calculate scale data")
    #   }
    # } else {
    sce <- Seurat::FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
    # }
  } else {
    if ("data" %in% layersNames) {
      print("load features from h5 file")
      Seurat::VariableFeatures(sce) <- rownames(h5_to_DF(h5[["var"]][["var"]]))
    } else {
      print("there is no features in h5 files and calFeatures is FALSE, please set 'calFeatures = TRUE' to calculate features")
    }
  }

  print("Add Scale...")
  # 如果calScale为TRUE，则计算scale数据，不读取data
  # 如果calScale为FALSE，且data存在，则读取data，不计算scale数据
  # 如果calScale为FALSE，且data存在，则读取data，判断dim是否一致，不一致则抛出提醒，是否calScale
  # 如果calScale为FALSE，且data不存在，则抛出提醒，是否calScale
  # Scale
  if (calScale) {
    if (calData) {
      sce <- Seurat::ScaleData(sce)
    } else {
      if ("data" %in% layersNames) {
        print("please set 'calData = TRUE' if you want to calculate scale data, load scale data from h5 file")
        scaleX <- h5_to_X(h5, assay, layer = "data", cellNames = cellNames, geneNames = geneNames)
        sce@assays[[assay]]@layers$scale.data <- NULL
        sce@assays[[assay]]@layers$scale.data <- scaleX
      } else {
        stop("there is no scale data in h5 files and calScale is FALSE,
          please set both 'calData = TRUE' and 'calScale = TRUE'
          to calculate scale data")
      }
    }
  } else {
    if ("data" %in% layersNames) {
      scaleX <- h5_to_X(h5, assay, "data")
      if (SeuratVersion == 5) {
        sce@assays[[assay]]@layers$scale.data <- NULL
        sce@assays[[assay]]@layers$scale.data <- scaleX
      } else {
        sce@assays[[assay]]@scale.data <- scaleX
      }

      if (all(dim(scaleX) != dim(sce))) {
        print("the dim of scale data is not equal to the dim of raw data, please set 'calScale = TRUE' to calculate new scale data")
      }
    } else {
      print("there is no scale data in h5 files and calScale is FALSE, please set both 'calData = TRUE' and 'calScale = TRUE' to calculate scale data")
    }
  }

  print("Add obs...")
  # obs
  sce <- Seurat::AddMetaData(sce, h5_to_obs(h5, "obs"))
  # graphs
  if ("graphs" %in% names(h5)) {
    sce <- sce_add_h5_to_graphs(sce, h5, cellNames, "graphs")
  }

  print("Add reductions...")
  # reductions
  if ("reductions" %in% names(h5)) {
    sce <- sce_add_h5_to_reductions(sce, h5, cellNames, assay, reductionsName = "reductions")
  }

  # # JoinLayers
  # if (SeuratVersion == 5) {
  #   if (!exists("scaleX") && all(dim(scaleX) == dim(sce))) {
  #     sce <- SeuratObject::JoinLayers(sce, assay = assay)
  #   }
  # }

  print("Add uns...")
  # 添加uns
  if ("uns" %in% names(h5)) {
    uns <- h5_to_uns(h5[["uns"]])
    for (unsName in names(uns)) {
      # class(uns[[unsName]]) <- "SeuratCommand"
      # attr(uns[[unsName]], "package") <- "SeuratObject"
      sce@commands[[unsName]] <- uns[[unsName]]
    }
  }

  print("Add images...")
  if ("images" %in% names(h5)) {
    images <- h5_to_images(h5[["images"]], assay, image_name = image_name, cellNames = cellNames)
    sce@images[[image_name]] <- images
  }
  return(sce)
}

h5_to_images <- function(images, assay, image_name = NULL, cellNames = NULL) {
  library(Seurat)
  coords <- as.data.frame(images[["coords"]][, ])
  if (!is.null(cellNames)) {
    rownames(coords) <- cellNames
  }
  colnames(coords) <- c("imagerow", "imagecol")
  image <- images[["image"]][, , ]
  scale_factors <- images[["scale_factors"]]
  fiducial <- scale_factors[["fiducial"]][]
  hires <- scale_factors[["hires"]][]
  lowres <- scale_factors[["lowres"]][]
  spot <- scale_factors[["spot"]][]

  scale.factors <- scalefactors(
    spot = spot,
    fiducial = fiducial,
    hires = hires,
    lowres = lowres
  )
  fov <- CreateFOV(as.data.frame(coords),
    type = "centroids", radius = spot,
    assay = assay, key = Key(image_name, quiet = TRUE)
  )
  visium.fov <- new(
    Class = "VisiumV2", boundaries = fov@boundaries,
    molecules = fov@molecules, assay = fov@assay, key = fov@key,
    image = image, scale.factors = scale.factors
  )
  # visium.fov@project.name <- image_name
  return(visium.fov)
}