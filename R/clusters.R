#script for functions to deal with centroids and cluster coords

# All clusters of circles are "clusterlist" ADTs:
# [["x"]] numeric vector of x coordinates of all circles
# [["y"]] numeric vector of y coordinates of all circles
# [["rad"]] numeric vector of radii of all circles
# [["centroid"]] numeric vector of the cluster centroid x and y coordinate
# [["clRad"]] approximate radius of the cluster
# [["clonotype"]] clonotype based on original clonecall

# getters for a single clusterlist
# does not handle the case of empty clusterlists

makeGetter <- function(key) {
    function(obj) if (isnt_empty(obj)) obj[[key]]
}

get_x_coords <- makeGetter(1)
get_y_coords <- makeGetter(2)
get_radii <- makeGetter(3)
get_centroid <- makeGetter(4)
get_centroid_x <- function(l) get_centroid(l)[1]
get_centroid_y <- function(l) get_centroid(l)[2]
get_cluster_radius <- makeGetter(5)
get_clonotypes <- function(l) l$clonotype # not index due to legacy clusterlists

get_num_clones <- function(l) length(get_x_coords(l))

contains_clonotypes <- function(x) !is.null(get_clonotypes(x))

# setters for a single clusterlist

makeSetter <- function(key) {
    function(obj, value) {
        if (isnt_empty(obj)) obj[[key]] <- value
        obj
    }
}

set_x_coords <- makeSetter(1)
set_y_coords <- makeSetter(2)
set_radii <- makeSetter(3)
set_centroid <- makeSetter(4)
set_cluster_radius <- makeSetter(5)
set_clonotypes <- makeSetter("clonotype")

# convert clusterlist to dataframe, assuming its ***valid***
# returns a dataframe with columns:
# - label
# - x
# - y
# - r
# - clonotype
convert_to_dataframe <- function(
    clstr_list, seurat_cluster_index, detail = TRUE
) {
    if (!detail) {
      return(data.frame(
        "label" = paste("cluster", seurat_cluster_index),
        "x" = get_centroid_x(clstr_list),
        "y" = get_centroid_y(clstr_list),
        "r" = get_cluster_radius(clstr_list)
      ))
    }

    data.frame(
        "label" = rep(
            paste("cluster", seurat_cluster_index),
            get_num_clones(clstr_list)
        ),
        "x" = get_x_coords(clstr_list),
        "y" = get_y_coords(clstr_list),
        "r" = get_radii(clstr_list)
    ) %>%
        update_clusterlist_df(clstr_list)
}

update_clusterlist_df <- function(df, clusterlist) {

    if (contains_clonotypes(df) || !contains_clonotypes(clusterlist)) return(df)
   
    if (nrow(df) == 0) {
        df$clonotype <- character(0)
        return(df)
    }

    df %>% dplyr::mutate(clonotype = get_clonotypes(clusterlist))

}

# centroid finder for a matrix of [x, y, cluster]
find_centroids <- function(xyc_df, ident_levels) {

  cll <- split(xyc_df, factor(xyc_df[, 3])) %>%
    lapply(function(x) c(mean(x[, 1]), mean(x[, 2])))
   
  list_output <- init_list(ident_levels, list())

  for (ident_level in ident_levels) {
    if (is.null(cll[[ident_level]])) next
    list_output[[ident_level]] <- cll[[ident_level]]
  }

  unname(list_output)
}

# get reduction centroids from seurat obj, where the barcodes in the reduction
# cell embeddings will be filtered to be exactly the same as those left in the
# metadata incase it was additionally filtered.
get_cluster_centroids <- function(seurat_obj, reduction, ident_levels) {

  if (!is_a_character(reduction)) {
    reduc_coords <- reduction@cell.embeddings[, 1:2]
  } else {
    reduc_coords <- get_2d_embedding(seurat_obj, reduction)
  }

  find_centroids(
    xyc_df = data.frame(
      rcppFilterReductionCoords(
        seuratBarcodes = rownames(seurat_obj@meta.data),
        reductionCoords = reduc_coords
      ),
      "ident" = seurat_obj@meta.data[["seurat_clusters"]]
    ),
    ident_levels
  )
}

# if new coord is null, assumes current coords are from the
# cluster being centered at (0, 0), and transforms x and y by
# centroid. Otherwise, translates by new_coord.
trans_coord <- function(cluster, new_coord = NULL) {

  if (!is.null(new_coord)) {
    trans_vec <- new_coord
    cluster <- set_centroid(cluster, get_centroid(cluster) + new_coord)
  } else {
    trans_vec <- get_centroid(cluster)
  }

  cluster %>%
    set_x_coords(get_x_coords(cluster) + trans_vec[1]) %>%
    set_y_coords(get_y_coords(cluster) + trans_vec[2])

}

# MOVE clusterlist to a new centroid, irrespective of previous centroid
move_cluster <- function(cluster, new_coord) {
  dx <- cluster[[4]][1] - new_coord[1]
  dy <- cluster[[4]][2] - new_coord[2]

  cluster[[1]] <- cluster[[1]] - dx
  cluster[[2]] <- cluster[[2]] - dy
  cluster[[4]] <- new_coord
  cluster
}

# function to GET a list of centroids from a list of clusterlists,
read_centroids <- function(list_of_clusterlists) {
  lapply(list_of_clusterlists, function(x) {
    if (is_empty(x)) list() else get_centroid(x)
  })
}

# clusterlist checkers below

isValidClusterList <- function(x, legacy = FALSE) {

    if (identical(x, list())) return(TRUE)

    basicChecks <- validate_that(
        is.list(x),
        (length(x) >= (6 - legacy)),
        is.numeric(x[[1]]),
        is.numeric(x[[2]]),
        is.numeric(x[[3]]),
        is.numeric(x[[4]]),
        is.numeric(x[[5]]),
        (legacy || is.character(x[[6]]))
    )

    if (!isTRUE(basicChecks)) {message(basicChecks); return(FALSE)}

    lengthChecks <- validate_that(
        length(unique(sapply(x[c(1, 2, 3, if (legacy) NULL else 6)], length))) == 1,
        length(x[[4]]) == 2,
        length(x[[5]]) == 1
    )
    
    if (!isTRUE(lengthChecks)) {message(lengthChecks); return(FALSE)}

    nameChecks <- identical(
        names(x)[1:5],
        c("x", "y", "rad", "centroid", "clRad")
    )
    if (!legacy) {
        nameChecks <- nameChecks && identical(names(x)[6], "clonotype")
    }

    if (isFALSE(nameChecks)) {
        message("names don't match: ", paste(names(x), collapse = ", "))
    }
    return(nameChecks)
}

isValidListOfClusterLists <- function(x, legacy = FALSE, verbose = FALSE) {
    for (i in seq_along(x)) {
        if (isValidClusterList(x[[i]], legacy)) next
        if (verbose) {
            message("clusterlist ", i, " is invalid")
        }
        return(FALSE)
    }
    return(TRUE)
}
