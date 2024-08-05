checkSeuratVersion <- function() {
    return(as.numeric(strsplit(as.character(packageVersion("Seurat")), split = "\\.")[[1]][1]))
}

# load
addDF <- function(oldDF, newDF, type = "row", reverse = FALSE) {
    if (type == "row") {
        if (nrow(oldDF) == 0) {
            oldDF <- newDF
        } else {
            if (reverse) {
                oldDF <- rbind(newDF, oldDF)
            } else {
                oldDF <- rbind(oldDF, newDF)
            }
        }
    } else if (type == "col") {
        if (ncol(oldDF) == 0) {
            oldDF <- newDF
        } else {
            if (reverse) {
                oldDF <- cbind(newDF, oldDF)
            } else {
                oldDF <- cbind(oldDF, newDF)
            }
        }
    } else {
        stop("Unknown type")
    }
    return(oldDF)
}

getShape <- function(h5data) {
    if ("shape" %in% hdf5r::h5attr_names(h5data)) {
        return(hdf5r::h5attr(h5data, "shape"))
    } else {
        return(calMatrixDim(h5data[["indptr"]][], h5data[["indices"]][]))
    }
}

getEncodingType <- function(h5data) {
    if ("encoding-type" %in% hdf5r::h5attr_names(h5data)) {
        return(hdf5r::h5attr(h5data, "encoding-type"))
    } else {
        return("unknown")
    }
}

calMatrixDim <- function(p, i) {
    return(c(length(p) - 1, max(i) + 1))
}

# save
h5AddAttribute <- function(h5data, name, value) {
    dtype <- NULL
    space <- NULL
    if (length(value) == 1) {
        space <- hdf5r::H5S$new("scalar")
        if (is.character(value)) {
            dtype <- hdf5r::H5T_STRING$new(type = "c", size = Inf)
            dtype$set_cset("UTF-8")
        } else if (is.logical(value)) {
            dtype <- hdf5r::h5types$H5T_NATIVE_HBOOL
            # dtype <- H5T$copy("H5T_NATIVE_HBOOL")
        }
    }

    h5data$create_attr(name, value, dtype = dtype, space = space)
}

h5dfAddIndex <- function(h5data, df, indexName = "_index") {
    h5data[[indexName]] <- rownames(df)
    h5AddAttribute(h5data[[indexName]], "encoding-type", "string-array")
    h5AddAttribute(h5data[[indexName]], "encoding-version", "0.2.0")
    # return(h5data)
}