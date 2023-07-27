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

# should probably fool proof this and the next function more
relocate_clusters <- function(seurat_obj, relocate_cluster, relocation_coord) {
  for(i in 1:length(relocate_cluster)) {
    cl_ind <- relocate_cluster[i] - 1
    seurat_obj@reductions[['apotc']]@clusters[[cl_ind]] <- move_cluster(
      seurat_obj@reductions[['apotc']]@clusters[[cl_ind]], relocation_coord[[i]]
    )
  }
  seurat_obj
}

nudge_clusters <- function(seurat_obj, nudge_cluster, nudge_vector) {
  for(i in 1:length(nudge_cluster)) {
    cl_ind <- nudge_cluster[i] - 1
    seurat_obj@reductions[['apotc']]@clusters[[cl_ind]] <- trans_coord(
      seurat_obj@reductions[['apotc']]@clusters[[cl_ind]], nudge_vector[[i]]
    )
  }
  seurat_obj
}

recolor_clusters <- function(seurat_obj, recolor_cluster, new_color) {
  for (i in 1:length(recolor_cluster)) {
    seurat_obj@reductions[['apotc']]@cluster_colors[recolor_cluster[i] - 1] <- new_color[[i]]
  }
  seurat_obj
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
# also possible to boot up a shiny window in the future?
AdjustAPOTC <- function(
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
  nudge_vector = NULL,

  recolor_cluster = -1L, #same as above
  new_color = NULL
) {

  # check for argument issues
  stop_str <- adjustAPOTC_stop_str(
    seurat_obj, new_rad_scale_factor, adjust_repulsion, max_repulsion_iter,
    relocate_cluster, relocation_coord, relocate_cluster, nudge_cluster,
    nudge_vector
  )
  if (!is.null(stop_str)) {
    stop(stop_str)
  }

  if (new_rad_scale_factor > 0) {
    seurat_obj <- change_rad_scale(seurat_obj, new_rad_scale_factor)
  }

  if (!identical(recolor_cluster, -1L)) {
    seurat_obj <- recolor_clusters(seurat_obj, recolor_cluster, new_color)
  }

  if (!identical(relocate_cluster, -1L)) {
    seurat_obj <- relocate_clusters(
      seurat_obj, relocate_cluster, relocation_coord
    )
  }

  if (!identical(nudge_cluster, -1L)) {
    seurat_obj <- nudge_clusters(seurat_obj, nudge_cluster, nudge_vector)
  }

  if (adjust_repulsion) {
    seurat_obj <- repulse_again(
      seurat_obj, repulsion_threshold, repulsion_strength, max_repulsion_iter,
      verbose
    )
  }

  seurat_obj
}

# it is probably super common to need to readjust the clone_scale_factor
change_clone_scale <- function(seurat_obj, new_clone_scale, verbose = TRUE) {

	if (is.null(seurat_obj@reductions[['apotc']])) {
		stop("please run RunAPOTC first")
	}

	if (verbose) {message("Repacking all clusters with new clone scale factor")}

	seurat_obj@reductions[['apotc']]@clone_scale_factor <- new_clone_scale

	seurat_obj@reductions[['apotc']]@clusters <- pack_into_clusterlists(
		sizes = get_processed_clone_sizes(seurat_obj@reductions[['apotc']]),
		centroids = seurat_obj@reductions[['apotc']]@centroids,
		num_clusters = seurat_obj@reductions[['apotc']]@num_clusters,
		rad_scale = seurat_obj@reductions[['apotc']]@rad_scale_factor,
		ORDER = get_cmd(seurat_obj, "ORDER"),
		scramble = get_cmd(seurat_obj, "scramble"),
		try_place = get_cmd(seurat_obj, "try_place"),
		verbose = verbose
	)

	seurat_obj
}
