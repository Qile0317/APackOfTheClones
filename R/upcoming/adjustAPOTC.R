# alter the radii of each circle in the seurat object
change_rad_scale <- function(seurat_obj, new_factor) {
  old_factor <- seurat_obj@reductions[['apotc']]@rad_scale_factor
  trans_factor <- new_factor / old_factor
  seurat_obj@reductions[['apotc']]@rad_scale_factor <- new_factor
  
  for (i in 1:seurat_obj@reductions[['apotc']]@num_clusters) {
    curr <- seurat_obj@reductions[['apotc']]@clusters[[i]]
    if (!isnt_empty(curr)) {
      next
    }
    seurat_obj@reductions[['apotc']]@clusters[[i]][[3]] <- curr[[3]] * trans_factor 
  }
  seurat_obj
}

# function for more concisce syntax for adjusting repulsion
repulse_again <- function(
  seurat_obj, repulsion_threshold, repulsion_strength, max_repulsion_iter, verbose
) {
  results <- get_repulsed_clusterlists_and_centroids(
    seurat_obj@reductions[['apotc']]@clusters,
    seurat_obj@reductions[['apotc']]@centroids,
    seurat_obj@reductions[['apotc']]@num_clusters,
    
    repulsion_threshold, repulsion_strength, max_repulsion_iter, verbose
  )
  seurat_obj@reductions[['apotc']]@clusters <- results[[1]]
  seurat_obj@reductions[['apotc']]@centroids <- results[[2]]
  seurat_obj
}

relocate_cluster <- function(seurat_obj, relocate_cluster, relocation_coord) {
  if (length(relocate_cluster) == length(relocation_coord)) {
    
  }
}

adjustAPOTC_stop_str <- function(
  seurat_obj, new_rad_scale_factor, adjust_repulsion, max_repulsion_iter,
  relocation_cluster, relocation_coord, relocate_cluster, nudge_cluster,
  nudge_vector
) {
  if (is.null(seurat_obj@reductions[['apotc']])) {
    return("RunAPOTC has not been ran on the input seurat object. Please follow the APackOfTheClones documentation and try again")
  }
  if (new_rad_scale_factor < 0) {
    return("new_rad_scale_factor must be a positive number")
  }
  if ((max_repulsion_iter < 1 || !is_int(max_repulsion_iter)) && adjust_repulsion) {
    return("max_repulsion_iter must be a positive integer")
  }
  if (!is.null(relocation_coord) && (length(relocate_cluster) != length(relocation_coord))) {
    return("length of relocate_cluster must be the same as the length of relocation_coord")
  }
  if (!is.null(nudge_vector) && (length(nudge_cluster) != length(nudge_vector))) {
    return("length of nudge_cluster must be the same as the length of nudge_vector")
  }
  if ((relocate_cluster != -1) && (nudge_cluster != -1)) {
    if (has_repeats(relocation_cluster, nudge_cluster)) {
      return("There are repeated elements in relocate_cluster and/or nudge_cluster")
    }
  }
  return(NULL)
}

# need functions for readjusting the apotc reduction for better visuals
# also possible to boot up a shiny window
adjustAPOTC <- function(
  seurat_obj,
  verbose = TRUE,
  
  new_rad_scale_factor = 0L, # assumes you dont want to change 
  
  adjust_repulsion = FALSE,
  repulsion_threshold = 1,
  repulsion_strength = 1,
  max_repulsion_iter = 10L,
  
  relocate_cluster = -1L, # can also be a vector
  relocation_coord = NULL, # vector or list of vectors
  
  nudge_cluster = -1L, # same as above
  nudge_vector = NULL
) {
  stop_str <- adjustAPOTC_stop_str(
    seurat_obj, new_rad_scale_factor, adjust_repulsion, max_repulsion_iter,
    relocation_cluster, relocation_coord, relocate_cluster, nudge_cluster,
    nudge_vector
  )
  if (!is.null(stop_str)) {
    stop(stop_str)
  }
  
  # need warnings first
  # # check if the seurat command is there
  # check if it was run correctly
  
  if (new_rad_scale_factor > 0) {
    seurat_obj <- change_rad_scale(seurat_obj, new_rad_scale_factor)
  }
  
  if (adjust_repulsion) {
    seurat_obj <- repulse_again(
      seurat_obj, repulsion_threshold, repulsion_strength, max_repulsion_iter,
      verbose
    )
  }
  
  seurat_obj
}