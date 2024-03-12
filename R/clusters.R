#script for functions to deal with centroids and cluster coords

# All clusters of circles are referred to as a "clusterlist" in code comments
# A clusterlist just refers to a simple R list of length five with:
# [["x"]] numeric vector of x coordinates of all circles
# [["y"]] numeric vector of y coordinates of all circles
# [["rad"]] numeric vector of radii of all circles
# [["centroid"]] numeric vector of the cluster centroid x and y coordinate
# [["clRad"]] approximate radius of the cluster
# [["clonotype"]] clonotype based on original clonecall

# centroid finder for a matrix of [x, y, cluster]
find_centroids <- function(df, total_clusters) {
  cll <- split(df, factor(df[, 3])) #the last cluster column becomes redundant
  l <- length(cll)

  nameset <- rep(c(""), times = l)
  xset <- rep(c(0), times = l)
  yset <- xset

  for (i in 1:l){
    nameset[i] <- cll[[i]][,3][1]
    xset[i] <- sum(cll[[i]][, 1]) / length(cll[[i]][, 1])
    yset[i] <- sum(cll[[i]][, 2]) / length(cll[[i]][, 2])
  }
  
  list_output <- init_list(num_elements = total_clusters, init_val = list())
  for (i in 1:l) {
    list_output[[as.integer(nameset[i])]] <- c(xset[i], yset[i])
  }
  return(list_output)
}

# get reduction centroids from seurat obj, where the barcodes in the reduction
# cell embeddings will be filtered to be exactly the same as those left in the
# metadata incase it was additionally filtered.
get_cluster_centroids <- function(
  seurat_obj, reduction = "umap", passed_in_reduc_obj = FALSE
) {

  if (passed_in_reduc_obj) {
    reduc_coords <- reduction@cell.embeddings[, 1:2]
  } else {
    reduc_coords <- get_2d_embedding(seurat_obj, reduction)
  }

  find_centroids(
    df = data.frame(
      rcppFilterReductionCoords(
        seuratBarcodes = rownames(seurat_obj@meta.data),
        reductionCoords = reduc_coords
      ),
      seurat_obj@meta.data[["seurat_clusters"]]
    ),
    total_clusters = get_num_total_clusters(seurat_obj)
  )
}

# TRANSFORM coordinates of a clusterlist from c(0, 0) to its own new centroid,
# or MOVE to new coord from current
trans_coord <- function(cluster, new_coord = NULL) {
  if (!is.null(new_coord)) {
    dx <- new_coord[1]
    dy <- new_coord[2]
    cluster[[4]] <- cluster[[4]] + c(dx, dy)
  } else {
    dx <- cluster[[4]][1]
    dy <- cluster[[4]][2]
  }
  cluster[[1]] <- cluster[[1]] + dx
  cluster[[2]] <- cluster[[2]] + dy
  cluster
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
  output_centroids <- init_list(length(list_of_clusterlists), list())
  for (i in seq_along(list_of_clusterlists)) {
    curr <- list_of_clusterlists[[i]]
    if (is.null(get_centroid(curr)) || is_empty(get_centroid(curr))) {
      next
    }
    output_centroids[[i]] <- get_centroid(curr)
  }
  output_centroids
}

# getters for a single clusterlist

get_x_coords <- function(l) l[[1]]
get_y_coords <- function(l) l[[2]]
get_radii <- function(l) l[[3]]
get_centroid <- function(l) l$centroid
get_cluster_radius <- function(l) l[[5]]
get_num_clones <- function(l) length(get_x_coords(l))

get_clonotypes <- function(x) UseMethod("get_clonotypes")
get_clonotypes.list <- function(x) x[["clonotype"]]
get_clonotypes.data.frame <- function(x) x[["clonotype"]]

contains_clonotypes <- function(x) !is.null(get_clonotypes(x))

# convert clusterlist to dataframe, assuming its valid
convert_to_dataframe <- function(clstr_list, seurat_cluster_index) {
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
