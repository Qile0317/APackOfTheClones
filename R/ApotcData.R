#' @title
#' The ApotcData class
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' An S4 class for storing information about T/B cell clonal expansion to be
#' used by various [APackOfTheClones] functions. Instances of this object type
#' are stored by [RunAPOTC] in a seurat object's `@misc` slot under a list named
#' `"APackOfTheClones"`. This class is not meant to be directly instantiated,
#' accessed, nor modified by the user. Attributes should only be accessed by the
#' associated getters so its independent of implementation.
#'
#' @slot reduction_base character indicating the reduction the plotting was
#' based off of.
#' @slot clonecall character indicating the column name of the seurat object's
#' metadata that contains the clone information.
#' @slot metadata_filter_string character indicating the metadata filter string
#' used to subset the seurat object before running APOTC.
#' @slot clusters Clustered clones which is an R list of lists of length 5,
#' with each list of length 5, with the first 3 elements being numeric vectors
#' containing the x and y coordinates and radii of each clone in the cluster.
#' The fourth element is a numeric vector of length 2 indicating the centroid of
#' the cluster, and the fifth is the estimated cluster radius.
#' @slot centroids Centroids of the clusters which default to the same centroids
#' of the clusterlists, tho this also means clusters has to be iterated over
#' everytime the plotting function is called. but it wont be slow probably
#' @slot clone_sizes the original unscaled clone sizes for the samples.
#' @slot clone_scale_factor scale factor to multiply `clone_sizes` by when
#' running the clonal expansion plotting algorithms
#' @slot rad_scale_factor scale factor to multiply the radii in clusterlists by
#' after they have all been computed to increase spacing between circles. Might
#' also be better in the future to instead just have a number to subtract :/
#' @slot cluster_colors character vector indicating coloration of each cluster
#' @slot reduction_base character indicating the reduction the plotting was
#' based off of
#' @slot label_coords list of numeric vectors of length two indicating the (x,y)
#' coordinates of each label if plotted
#'
#' @keywords internal
#' @noRd
#'
methods::setClass(
	Class = "ApotcData",
	slots = c(
		reduction_base = 'character',
		clonecall = 'character',
		metadata_filter_string = 'character',

		clusters = 'list',
		centroids = 'list',
		clone_sizes = 'list',
		num_clusters = 'numeric',

		clone_scale_factor = 'numeric', # should make these overridable in APOTCPlot, and choose if they want to modify
		rad_scale_factor = 'numeric',
		cluster_colors = 'character',
		labels = 'character',
		label_coords = 'list'
	)
)

ApotcData <- function(
	seurat_obj, metadata_filter_condition, clonecall, reduction_base,
	clone_scale_factor, rad_scale_factor
) {
	if (identical(metadata_filter_condition, "")) {
		return(initializeApotcData(
			seurat_obj, clonecall, reduction_base, clone_scale_factor, rad_scale_factor
		))
	}
	initializeSubsetApotcData(
		seurat_obj, metadata_filter_condition, clonecall, reduction_base,
		clone_scale_factor, rad_scale_factor
	)
}

initializeApotcData <- function(
	seurat_obj, clonecall, reduction_base, clone_scale_factor, rad_scale_factor
) {
	num_clusters <- get_num_total_clusters(seurat_obj) # may need to redo, will leave empty elements for invalid clusters
	initial_centroids <- get_cluster_centroids(seurat_obj, reduction_base)
	raw_all_clone_sizes <- count_raw_clone_sizes(
		seurat_obj = seurat_obj, num_clusters = num_clusters, clonecall = clonecall
	)

	methods::new(
		Class = "ApotcData",

		reduction_base = reduction_base,
		clonecall = clonecall,
		metadata_filter_string = "",

		clusters = list(),
		centroids = initial_centroids,
		clone_sizes = raw_all_clone_sizes,
		num_clusters = num_clusters,

		clone_scale_factor = clone_scale_factor,
		rad_scale_factor = rad_scale_factor,
		cluster_colors = gg_color_hue(num_clusters),
		labels = gen_labels(num_clusters),
		label_coords = initial_centroids
	)
}

# create subset / new obj based on new conditions, assuming valid!
# based on an initialized (not nessecarily packed) apotc obj
# assumes metadata_filter_condition cannot be null

initializeSubsetApotcData <- function(
	seurat_obj, metadata_filter_condition, clonecall, reduction_base,
	clone_scale_factor, rad_scale_factor
) {
	# subset the seurat metadata
	seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::filter(eval(parse(
		text = metadata_filter_condition
	)))

	apotc_obj <- initializeApotcData(
		seurat_obj, clonecall, reduction_base, clone_scale_factor, rad_scale_factor
	)

	apotc_obj@metadata_filter_string <- metadata_filter_condition
	apotc_obj
}

# pack the clones assuming centroids are present
circlepackClones <- function(apotc_obj, ORDER, scramble, try_place, verbose) {

	apotc_obj@clusters <- pack_into_clusterlists(
		sizes = get_processed_clone_sizes(apotc_obj),
		centroids = get_centroids(apotc_obj),
		num_clusters = get_num_clusters(apotc_obj),
		rad_decrease = get_rad_decrease(apotc_obj),
		ORDER = ORDER,
		scramble = scramble,
		try_place = try_place,
		verbose = verbose
	)

	# see which elemens of sizes are empty and set corresponding elements empty
	for (i in seq_len(get_num_clusters(apotc_obj))) {
		if (isnt_empty(apotc_obj@clusters[[i]])) {
			next
		}
		apotc_obj@centroids[[i]] <- list()
		apotc_obj@label_coords[[i]] <- list()
		# technically colors too :/
	}

	apotc_obj
}

# function to do repulsion for both cases

repulseClusters <- function(
	apotc_obj, repulsion_threshold, repulsion_strength, max_repulsion_iter,
	verbose
) {
	repulsed_clusters <- get_repulsed_clusterlists(
	    packed_clusters = get_clusterlists(apotc_obj),
	    initial_centroids = get_centroids(apotc_obj),
	    num_clusters = get_num_clusters(apotc_obj),
	    repulsion_threshold = repulsion_threshold,
		repulsion_strength = repulsion_strength,
		max_repulsion_iter = max_repulsion_iter,
		verbose = verbose
	)

	setModifiedClusterlists(
		apotc_obj, modified_clusterlists = repulsed_clusters
	)
}

# function to modify the apotc_obj's relevant slots when modified clusterlists
# are introduced e.g. for cluster repulsion or relocation. This cannot be used
# for completely new irrelevant clusterlists, as the centroids and label_coords
# are modified correspondingly to the original clusters.
setModifiedClusterlists <- function(apotc_obj, modified_clusterlists) {

	original_centroids <- get_centroids(apotc_obj)
	modified_centroids <- read_centroids(modified_clusterlists)

	apotc_obj@label_coords <- operate_on_same_length_lists(
        func = add,
        l1 = get_label_coords(apotc_obj),
        l2 = operate_on_same_length_lists(
            func = subtract,
            l1 = modified_centroids,
            l2 = original_centroids
        )
    )

	apotc_obj@clusters <- modified_clusterlists
	apotc_obj@centroids <- modified_centroids

	apotc_obj
}

convert_to_rad_decrease <- function(clone_scale_factor, rad_scale_factor) {
	clone_scale_factor * (1 - rad_scale_factor)
}

# getters

get_reduction_base <- function(apotc_obj) {
	apotc_obj@reduction_base
}

get_clonecall <- function(apotc_obj) {
	apotc_obj@clonecall
}

get_metadata_filter_string <- function(apotc_obj) {
	apotc_obj@metadata_filter_string
}

get_clusterlists <- function(apotc_obj) {
	apotc_obj@clusters
}

get_centroids <- function(apotc_obj) {
	apotc_obj@centroids
}

get_raw_clone_sizes <- function(apotc_obj) {
	apotc_obj@clone_sizes
}

get_processed_clone_sizes <- function(apotc_obj) {
  raw_tabled_clone_sizes <- get_raw_clone_sizes(apotc_obj)
  processed_sizes <- init_list(get_num_clusters(apotc_obj), list())

  for (i in seq_len(get_num_clusters(apotc_obj))) {
    if (!is_empty_table(raw_tabled_clone_sizes[[i]])) {
      processed_sizes[[i]] <- apotc_obj@clone_scale_factor *
        sqrt(as.numeric(raw_tabled_clone_sizes[[i]][[1]]))
    }
  }
  processed_sizes
}

get_num_clones <- function(apotc_obj) { # TODO test
	sum(unlist(get_raw_clone_sizes(apotc_obj)))
}

get_num_clusters <- function(apotc_obj) {
	apotc_obj@num_clusters
}

get_valid_num_clusters <- function(apotc_obj) {
	n <- 0
	for (cluster in apotc_obj@clusters) {
		if (isnt_empty(cluster)) {
			n <- n + 1
		}
	}
	n
}

get_clone_scale_factor <- function(apotc_obj) {
	apotc_obj@clone_scale_factor
}

get_rad_scale_factor <- function(apotc_obj) {
	apotc_obj@rad_scale_factor
}

get_rad_decrease <- function(apotc_obj) {
	convert_to_rad_decrease(
		clone_scale_factor = get_clone_scale_factor(apotc_obj),
		rad_scale_factor = get_rad_scale_factor(apotc_obj)
	)
}

get_cluster_colors <- function(apotc_obj) {
	apotc_obj@cluster_colors
}

get_labels <- function(apotc_obj) {
	apotc_obj@labels
}

get_label_coords <- function(apotc_obj) {
	apotc_obj@label_coords
}
