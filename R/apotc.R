#' The apotc (APackOfTheClones) reduction class
#'
#' A class for storing information about T cell clonal expansion in the seurat
#' `reductions` attribute
#'
#' @slot clusters The original circle packing "reduction" with no modifications,
#' a list of clusterlists that may include NAs.
#' @slot centroids Centroids of the clusters which default to the same centroids
#' of the clusterlists, tho this also means clusters has to be iterated over
#' everytime the plotting function is called. but it wont be slow probably
#' @slot clone_sizes the original unscaled clone sizes. Not sure if I should
#' leave it named
#' @slot clone_scale_factor self explanatory
#' @slot rad_scale_factor self explanatory
#'
#' @noRd
#'
setClass(
  Class = "apotc",
  slots = c(
    clusters = 'list',
    num_clusters = 'numeric',
    centroids = 'list', 
    clone_sizes = 'list', 
    clone_scale_factor = 'numeric',
    rad_scale_factor = 'numeric'
  )
)

# impl in rust or c++ and see if theres speed improvement
get_centroid_list <- function(clusterlists, num_clusters) {
  centroid_list <- vector(mode = "list", length = num_clusters)
  for (i in 1:num_clusters) {
    if (!is.na(clusterlists[[i]])) {
      centroid_list[[i]] <- clusterlists[[i]][[5]]
    }else {
      centroid_list[[i]] <- NA
    }
  }
  centroid_list
}

create_apotc <- function(
    clusters,
    clone_sizes,
    centroids = NULL,
    clone_scale_factor = 1,
    rad_scale_factor = 1) {
  
  new(
    Class = 'apotc',
    clusters = clusters, 
    centroids = NULL, 
    clone_sizes = clone_sizes, 
    clone_scale_factor = clone_scale_factor,
    rad_scale_factor = rad_scale_factor
  )
}

# there needs to be a RunAPOTC() and tune_apotc_param() or something like that
# # unginished script