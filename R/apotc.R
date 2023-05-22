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
#' running the clonal expansion plotting algorithmd
#' @slot rad_scale_factor scale factor to multiply the radii in clusterlis  ts by
#' after they have all been computed to incr ease spacing between circles. Might
#' also be better in the future to instead just have a number to subtract :/
#' @slot cluster_colors character vector indicating coloration of each cluster
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
    cluster_colors = 'character' 
  )
)

# initialize the reduction object from precomputed clusterlists
initialize_apotc <- function(
    num_clusters, clone_scale_factor = 1, rad_scale_factor = 1
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
    cluster_colors = gg_color_hue(num_clusters)
  )
}

# run all the packing algorithms and modify the apotc object
# should be ran after initialize_apotc
run_packing_algos <- function(
  apotc_obj, integrated_seurat_obj, clone_scale_factor,
  rad_scale_factor = 1, ORDER = TRUE,
  try_place = FALSE, verbose = TRUE,
  repulse = FALSE,
  repulsion_threshold = 1,
  repulsion_strength = 1,
  max_repulsion_iter = 10
) {
  apotc_obj <- add_raw_clone_sizes(apotc_obj, integrated_seurat_obj)
  initial_centroids <- get_cluster_centroids(integrated_seurat_obj)
  
  packed_clusters <- pack_into_clusterlists(
    multiply_all(
      apotc_obj@clone_sizes, apotc_obj@num_clusters, clone_scale_factor
    ),
    initial_centroids,
    apotc_obj@num_clusters,
    rad_scale_factor,
    ORDER,
    try_place,
    verbose
  )
  
  if (repulse) {
    if(verbose){
      message(
        paste("\nrepulsing all clusters | max iterations =", max_repulsion_iter)
      )
    }
    packed_clusters <- repulse_cluster(
      packed_clusters, repulsion_threshold, repulsion_strength,
      max_repulsion_iter, verbose
    )
    initial_centroids <- read_centroids(
      packed_clusters, initial_centroids, apotc_obj@num_clusters
    )
  }
  
  apotc_obj@clusters <- packed_clusters
  apotc_obj@centroids <- initial_centroids
  
  apotc_obj
}

# alias to add the apotc object to the seurat obj for readability
combined_obj <- function(seurat_obj, apotc_obj) {
  seurat_obj@reductions[['apotc']] <- apotc_obj
  seurat_obj
}

# function to imitate RunUMAP
RunAPOTC <- function(
  seurat_obj,
  tcr_df = "seurat_obj_already_integrated",
  reduction_base = "umap",
  clone_scale_factor = 0.1, # do 0.5 for test ds - need to make an estimator based on testing
  rad_scale_factor = 0.95, 
  ORDER = TRUE,
  try_place = FALSE,
  verbose = TRUE,
  repulse = FALSE, 
  repulsion_threshold = 1,
  repulsion_strength = 1,
  max_repulsion_iter = 10) {
  
  # add seurat command
  seurat_obj@commands[["RunAPOTC"]] <- make_apotc_command()
  
  # errors/warnings:
  if (is.null(seurat_obj@reductions[[reduction_base]])) {
    stop(paste("No", reduction_base, "reduction found on the seurat object"))
  }
  if ((!is.data.frame(tcr_df)) && is.null(seurat_obj@meta.data[["raw_clonotype_id"]])) {
    stop("Seurat object is missing the raw_clonotype_id data or isn't integrated with the TCR library. Consider integrating the T-cell library into the seurat object again.")
  }
  if (max_repulsion_iter > 1000) {
    warning("Repulsion iteration count > 1000, consider reducing max_repulsion_iter if runtime is too long")
  }
  
  # integrate TCR - in future remake to be compatible with scRepertoire
  if (is.data.frame(tcr_df)) {
    seurat_obj <- integrate_tcr(seurat_obj, tcr_df, verbose = verbose)
  }
  
  # get number of seurat clusters
  num_clusters <- get_num_clusters(seurat_obj)
  
  # initialize apotc S4 class
  apotc_obj <- initialize_apotc(
    num_clusters, clone_scale_factor, rad_scale_factor
  )
  
  # run all packing algos 
  apotc_obj <- run_packing_algos(
    apotc_obj, seurat_obj, rad_scale_factor, ORDER, try_place, verbose,
    repulse, repulsion_threshold, repulsion_strength, max_repulsion_iter
  )

  # add the finished apotc object to reductions and return
  seurat_obj <- combined_obj(seurat_obj, apotc_obj)
  seurat_obj
}
