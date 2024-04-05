#script for functions to deal with centroids and cluster coords

# All clusters of circles are "clusterlist" ADTs:
# [["x"]] numeric vector of x coordinates of all circles
# [["y"]] numeric vector of y coordinates of all circles
# [["rad"]] numeric vector of radii of all circles
# [["centroid"]] numeric vector of the cluster centroid x and y coordinate
# [["clRad"]] approximate radius of the cluster
# [["clonotype"]] clonotype based on original clonecall

# getters for a single clusterlist

get_x_coords <- function(l) l[[1]]
get_y_coords <- function(l) l[[2]]
get_radii <- function(l) l[[3]]
get_centroid <- function(l) l[[4]]
get_cluster_radius <- function(l) l[[5]]
get_clonotypes <- function(l) l$clonotype # not index due to legacy clusterlists

get_num_clones <- function(l) length(get_x_coords(l))

contains_clonotypes <- function(x) !is.null(get_clonotypes(x))

# setters for a single clusterlist
set_x_coords <- function(l, v) {l$x <- v; l}
set_y_coords <- function(l, v) {l$y <- v; l}
set_radii <- function(l, v) {l$rad <- v; l}
set_centroid <- function(l, v) {l$centroid <- v; l}
set_cluster_radius <- function(l, v) {l$clRad <- v; l}
set_clonotypes <- function(l, v) {l$clonotype <- v; l}

# convert clusterlist to dataframe, assuming its valid
# TODO allow get_abstract option to just get 1 large circle
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
