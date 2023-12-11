#' @title
#' The ApotcData class (TO CHANGE)
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' An S4 class for storing information about T/B cell clonal expansion to be
#' used by various [APackOfTheClones] functions. Instances of this object type
#' are stored by [RunAPOTC] in a seurat object's `@misc` slot under a list named
#' `"APackOfTheClones"`. This class is not meant to be directly instantiated nor
#' modified by the user. Read the vignette about the various getters for this
#' class.
#'
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
#' @export
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

		clone_scale_factor = 'numeric', # make these overridable in APOTCPlot, and choose if they want to modify
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
# creates a temporary apotc obj if it doesn exist already
# assumes metadata_filter_condition cannot be null

initializeSubsetApotcData <- function(
	seurat_obj, metadata_filter_condition, clonecall, reduction_base,
	clone_scale_factor, rad_scale_factor
) {
	# subset the seurat metadata
	seurat_obj@meta.data %>% dplyr::filter(eval(parse(
		text = metadata_filter_condition
	)))

	apotc_obj <- initializeApotcData(
		seurat_obj, clonecall, reduction_base, clone_scale_factor, rad_scale_factor
	)

	apotc_obj@metadata_filter_string <- metadata_filter_condition
	apotc_obj
}

circlepackClones <- function(apotc_obj, ORDER, scramble, try_place, verbose) {

	rad_decrease <- convert_to_rad_decrease(
		apotc_obj@rad_scale_factor, apotc_obj@clone_scale_factor
	)

	apotc_obj@clusters <- pack_into_clusterlists(
		sizes = get_processed_clone_sizes(apotc_obj),
		centroids = apotc_obj@centroids,
		num_clusters = apotc_obj@num_clusters,
		rad_decrease = rad_decrease,
		ORDER = ORDER,
		scramble = scramble,
		try_place = try_place,
		verbose = verbose
	)

	apotc_obj
}

# function to do repulsion for both cases

repulseClusters <- function(
	apotc_obj, repulsion_threshold, repulsion_strength, max_repulsion_iter,
	verbose
) {
	repulsion_results <- get_repulsed_clusterlists_and_centroids(
		apotc_obj@clusters, apotc_obj@centroids,
		repulsion_threshold, repulsion_strength, max_repulsion_iter, verbose
	)

	apotc_obj@clusters <- repulsion_results[[1]]
	apotc_obj@centroids <- repulsion_results[[2]]
	apotc_obj@label_coords <- repulsion_results[[2]] # TODO if label coords were modded perhaps they should only move by a factor instead?

	apotc_obj
}



# should have getters and setters
#
# # should have a matchcolors func that takes into accoutn cluster names
