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
		sample_prefixes = 'character',

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

# internal constructor for the .defaultApotcDataSample sample case

ApotcData <- function(
	seurat_obj, clone_scale_factor, rad_scale_factor, reduction_base
) {
	empty_list <- vector("list", num_clusters)
	methods::new(
		Class = 'ApotcData',

		reduction_base = reduction_base,
		clusters = empty_list,
		centroids = empty_list,
		clone_sizes = empty_list,
		num_clusters = num_clusters,

		clone_scale_factor = clone_scale_factor,
		rad_scale_factor = rad_scale_factor,
		cluster_colors = gg_color_hue(num_clusters),
		labels = gen_labels(num_clusters),
		label_coords = empty_list
	)
}

.defaultApotcDataSample <- "__all__"
utils::globalVariables(c(".defaultApotcDataSample"))

# new format, there will be a list of apotc objects in the seurat@misc slot. the list will be named apotc.
# each is dependent on reduction/samples and within the list there will be named elements for each reduction/sample combo
# and make it optional in RunAPOTC if this should be stored. APOTCPlot will then be able to have the apotc obj slot input
# alternatively the sample/id configuration.

# should have getters and setters
#
# # should have a matchcolors func that takes into accoutn cluster names
# should by default force the "all" run for the first time even if it isnt selected



# all subsequent apotc data that isnt "all" should be derived from
# the "all" sample case with the corresponding reduction.

# warn helper
run_apotc_warn_str <- function(
		seurat_obj, reduction_base, clone_scale_factor, ORDER, scramble
) {
	if (tolower(reduction_base) == 'apotc') {
		return("please only use the umap, tsne, or pca reduction")
	}
	if (is.null(seurat_obj@reductions[[attempt_correction(reduction_base)]])) {
		return(paste(
			"No", reduction_base, "reduction found on the seurat object,",
			"ensure the the reduction has been computed. Otherwise, did you",
			"mean", closest_word(reduction_base), "?"
		))
	}
	if (should_estimate(clone_scale_factor) && (clone_scale_factor <= 0)) {
		return("clone_scale_factor has to be a positive number")
	}
	if (ORDER && scramble) {
		return("ORDER and scramble are both TRUE, please set only one to TRUE")
	}
	return(NULL)
}
