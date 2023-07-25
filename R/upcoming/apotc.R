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
#' @slot clone_scale_factor scale factor to multiply `clone_sizes` by when
#' running the clonal expansion plotting algorithms
#' @slot rad_scale_factor scale factor to multiply the radii in clusterlists by
#' after they have all been computed to increase spacing between circles. Might
#' also be better in the future to instead just have a number to subtract :/
#' @slot cluster_colors character vector indicating coloration of each cluster
#' @slot reduction_base character indicating the reduction the plotting was
#' based off of
#'
#' @keywords internal
#'
setClass(
  Class = "apotc",
  slots = c(
    clusters = 'list',
    num_clusters = 'numeric',
    centroids = 'list',
    clone_sizes = 'list',
    clone_scale_factor = 'numeric',
    rad_scale_factor = 'numeric',
    cluster_colors = 'character',
    reduction_base = 'character'
  )
)

# initialize the reduction object from precomputed clusterlists
initialize_apotc <- function(
  num_clusters, clone_scale_factor, rad_scale_factor, reduction_base
) {
  empty_list <- vector("list", num_clusters)
  methods::new(
    Class = 'apotc',
    clusters = empty_list,
    num_clusters = num_clusters,
    centroids = empty_list,
    clone_sizes = empty_list,
    clone_scale_factor = clone_scale_factor,
    rad_scale_factor = rad_scale_factor,
    cluster_colors = gg_color_hue(num_clusters),
    reduction_base = reduction_base
  )
}

# function to imitate RunUMAP - something is wrong w the C++ packing in this function
# session aborts on cluster 2 during testing
RunAPOTC <- function(
  seurat_obj,
  tcr_df = "seurat_obj_already_integrated",
  reduction_base = "umap",
  clone_scale_factor = 0.1, # do 0.5 for test ds - need to make an estimator based on testing
  rad_scale_factor = 0.95,
  ORDER = TRUE,
  scramble = FALSE,
  try_place = FALSE,
  verbose = TRUE,
  repulse = FALSE,
  repulsion_threshold = 1,
  repulsion_strength = 1,
  max_repulsion_iter = 10L
) {
  # call time for seurat commands
  call_time = Sys.time()

  reduction_base <- attempt_correction(reduction_base)

  # errors/warnings:
  if (is.null(seurat_obj@reductions[[reduction_base]])) {
    stop(paste(
        "No", reduction_base, "reduction found on the seurat object,",
        "please ensure spelling is correct"
    ))
  }
  if ((!is.data.frame(tcr_df)) && is.null(seurat_obj@meta.data[["raw_clonotype_id"]])) {
    stop("Seurat object is missing the raw_clonotype_id data or isn't integrated with the TCR library. Consider integrating the T-cell library into the seurat object again.")
  }
  if (clone_scale_factor <= 0) {
    stop("clone_scale_factor has to be a positive number")
  }

  if (verbose) {message("Initializing APOTC run")}

  # integrate TCR - in future remake to be compatible with scRepertoire
  if (is.data.frame(tcr_df)) {
      seurat_obj <- dev_integrate_tcr(seurat_obj, tcr_df, verbose, call_time)
  }

  # add seurat command
  seurat_obj@commands[["RunAPOTC"]] <- make_apotc_command(call_time)

  # initialize apotc S4 class
  apotc_obj <- initialize_apotc(
    num_clusters = get_num_clusters(seurat_obj),
    clone_scale_factor,
    rad_scale_factor,
    reduction_base
  )

  # infer clone sizes and centroids
  apotc_obj <- add_raw_clone_sizes(apotc_obj, seurat_obj)
  initial_centroids <- get_cluster_centroids(
    seurat_obj, reduction = reduction_base
  )

  # pack the clusterlists
  packed_clusters <- pack_into_clusterlists(
    sizes = get_processed_clone_sizes(apotc_obj),
    centroids = initial_centroids,
    num_clusters = apotc_obj@num_clusters,
    rad_scale = rad_scale_factor,
    ORDER = ORDER,
    scramble = scramble,
    try_place = try_place,
    verbose = verbose
  )

  if (repulse) {
    results <- get_repulsed_clusterlists_and_centroids(
      packed_clusters, initial_centroids, apotc_obj@num_clusters,
      repulsion_threshold, repulsion_strength, max_repulsion_iter, verbose
    )
    packed_clusters <- results[[1]]
    initial_centroids <- results[[2]]
  }

  # add the clusterlists and centroids to the obj
  apotc_obj@clusters <- packed_clusters
  apotc_obj@centroids <- initial_centroids

  # add the finished apotc object to reductions, print message, and return
  seurat_obj@reductions[['apotc']] <- apotc_obj
  if (verbose) {print_completion_time(call_time)}
  seurat_obj
}
