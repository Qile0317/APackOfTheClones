#' @title
#' Adjust the parameters of the APackOfTheClones reduction in a seurat
#' object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' If the user is unsatisfied with the clonal expansion plot that
#' was generated from `RunAPOTC` and `APOTCPlot`, this function has a range of
#' arguments to modify the data and/or parameters of the visualization. Note
#' that some of the arguments may conflict with eachother.
#'
#' @inheritParams RunAPOTC
#'
#' @param seurat_obj The seurat object to be adjusted.
#' @param new_clone_scale_factor a single numeric in (0, 1]. changes the
#' clone_scale_factor
#' @param new_rad_scale_factor a single numeric in (0, 1]. changes the
#' radius scaling factor of all circles.
#' @param relocate_cluster numeric of arbitrary length. Indicates which
#' cluster(s) to relocate to new coordinates
#' @param relocation_coord numeric of length two or a list of numerics of length
#' two of length of `relocate_cluster`. If its a list, indicates each coordinate
#' that the clusters in `relocate_cluster` should move to. If its just a
#' numeric, then will relocate all clusters in `relocate_cluster` to the input,
#' which is likely not desired behavior, so this should only be convenience
#' syntax if `relocate_cluster` has length 1.
#' @param nudge_cluster numeric of arbitrary length. Indicates which
#' cluster(s) to "nudge"/translate their coordinate(s) by.
#' @param nudge_vector numeric of length two or a list of numerics of length
#' two of length of `nudge_cluster`. If its a list, indicates each translation
#' vector (in other words, x-y coordinates) that the clusters in
#' `nudge_cluster` should be translate by. If its just a numeric, then will
#' translate all clusters in `nudge_cluster` by the input - which mostly is
#' syntactic sugar for translating a single cluster if the input of
#' `nudge_cluster` is of length 1.
#' @param recolor_cluster numeric of arbitrary length. Indicates which
#' cluster(s) to change their color by.
#' @param new_color character of arbitrary length. Indicates the corresponding
#' new colors that selected clusters in `recolor_cluster` should be changed to.
#' @param rename_label Numeric or character. Indicates the index or name of
#' label(s) to be renamed.
#' @param new_label Character. Indicates the corresponding new label(s) that
#' selected label(s) in `rename_label` should be changed to.
#' @param relocate_label Numeric or character. Indicates the index or name of
#' label(s) to be relocated.
#' @param label_relocation_coord Numeric of length two or a list of numerics of
#' length two of length of `relocate_label`. If it's a list, indicates each
#' coordinate that the labels in `relocate_label` should move to. If it's just
#' a numeric, then will relocate all labels in `relocate_label` to the input,
#' which is likely not desired behavior, so this should only be convenience
#' syntax if `relocate_label` has length 1.
#' @param nudge_label Numeric or character. Indicates the index or name of
#' label(s) to be "nudged"/translated.
#' @param label_nudge_vector Numeric of length two or a list of numerics of
#' length two of length of `nudge_label`. If it's a list, indicates each
#' translation vector (in other words, x-y coordinates) that the labels in
#' `nudge_label` should be translated by. If it's just a numeric, then will
#' translate all labels in `nudge_label` by the input - which mostly is
#' syntactic sugar for translating a single label if the input of
#' `nudge_label` is of length 1.
#'
#' @return The adjusted `seurat_obj`
#' @export
#'
#' @examples
#' # do an APackOfTheClones run
#' pbmc <- RunAPOTC(get(data("combined_pbmc")), verbose = FALSE)
#'
#' # adjust the rad_scale_factor, and nudge cluster 1 by (1, 1)
#' pbmc <- AdjustAPOTC(
#'     pbmc,
#'     new_rad_scale_factor = 0.9,
#'     nudge_cluster = 1,
#'     nudge_vector = c(1, 1),
#'     verbose = FALSE
#' )
#'
#' # plot the result
#' APOTCPlot(pbmc)
#'
#' # perhaps multiple clusters need to be relocated and repulsed
#' pbmc <- AdjustAPOTC(
#'     pbmc,
#'     relocate_cluster = c(5, 10),
#'     relocation_coord = list(c(2, 3.5), c(0, 5)),
#'     repulse = TRUE,
#'     verbose = FALSE
#' )
#'
#' # plot again to check results
#' APOTCPlot(pbmc)
#'
AdjustAPOTC <- function(
	seurat_obj,
	reduction_base = NULL,
	clonecall = NULL,
	...,
	extra_filter = NULL,
	run_id = NULL,

	new_rad_scale_factor = NULL,
	new_clone_scale_factor = NULL,

	repulse = FALSE,
	repulsion_threshold = 1,
	repulsion_strength = 1,
	max_repulsion_iter = 10L,

	relocate_cluster = NULL,
	relocation_coord = NULL,
	nudge_cluster = NULL,
	nudge_vector = NULL,
	recolor_cluster = NULL,
	new_color = NULL,

	rename_label = NULL, # TODO add functionality to change all at once or modify all at once
	new_label = NULL,
	relocate_label = NULL,
	label_relocation_coord = NULL,
	nudge_label = NULL,
	label_nudge_vector = NULL,

	# # upcoming
	# interactive = FALSE,
	verbose = TRUE
) {
	varargs_list <- list(...)
	args <- hash::hash(as.list(environment()))
	AdjustAPOTC_error_handler(args)

	args$object_id <- infer_object_id_if_needed(args, varargs_list = varargs_list)
	apotc_obj <- getApotcData(seurat_obj, args$object_id)

	# # TODO
	# if (interactive) {
	# 	return(runShinyAdjustAPOTC(args))
	# }

	if (should_change(new_clone_scale_factor)) {
		apotc_obj <- change_clone_scale(seurat_obj, args)
	}

	if (should_change(new_rad_scale_factor)) {
		apotc_obj <- change_rad_scale(apotc_obj, new_rad_scale_factor)
	}

	if (should_change(recolor_cluster)) {
		apotc_obj <- recolor_clusters(apotc_obj, recolor_cluster, new_color)
	}

	if (should_change(relocate_cluster)) {
		apotc_obj <- relocate_clusters(
			apotc_obj, relocate_cluster, relocation_coord
		)
	}

	if (should_change(nudge_cluster)) {
		apotc_obj <- nudge_clusters(apotc_obj, nudge_cluster, nudge_vector)
	}

	if (repulse) {
		apotc_obj <- repulseClusters(
			apotc_obj, repulsion_threshold, repulsion_strength,
			max_repulsion_iter, verbose
		)
		if (verbose) message()
	}

	if (should_change(rename_label)) {
		apotc_obj <- rename_labels(apotc_obj, rename_label, new_label)
	}

	if (should_change(relocate_label)) {
		apotc_obj <- relocate_labels(
			apotc_obj, relocate_label, label_relocation_coord
		)
	}

	if (should_change(nudge_label)) {
		apotc_obj <- nudge_labels(apotc_obj, nudge_label, label_nudge_vector)
	}

	setApotcData(seurat_obj, args$object_id, apotc_obj)
}

	# rename_label = NULL, # TODO add functionality to change all at once or modify all at once
	# new_label = NULL,
	# relocate_label = NULL,
	# label_relocation_coord = NULL,
	# nudge_label = NULL,
	# label_nudge_vector = NULL,

AdjustAPOTC_error_handler <- function(args) {

	check_apotc_identifiers(args)

	args$new_rad_scale_factor %>% typecheck(is_a_positive_numeric, is.null)
	args$new_clone_scale_factor %>% typecheck(is_a_positive_numeric, is.null)

	check_repulsion_params(args)

	# TODO check lengths for a bunch of these

	args$relocate_cluster %>% typecheck(is_integer, is.null)
	args$relocation_coord %>% typecheck(
		is_numeric_pair, is_list_of_numeric_pair, is.null
	) 

	# if (!is.null(args$relocation_coord) && (length(args$relocate_cluster) != length(args$relocation_coord))) {
	# return("length of relocate_cluster must be the same as the length of relocation_coord")
	# }

	args$nudge_cluster %>% typecheck(is_integer, is.null)
	args$nudge_vector %>% typecheck(
		is_numeric_pair, is_list_of_numeric_pair, is.null
	)

	args$recolor_cluster %>% typecheck(is_integer, is.null)
	args$new_color %>% typecheck(is_character, is.null)

	# TODO more checks
}

# TODO it should be more efficient + safer to mathematically transform all vals
change_clone_scale <- function(seurat_obj, args) {

	if (args$verbose) {
		message("Repacking all clusters with new clone scale factor")
	}

	past_params <- find_seurat_command(
		seurat_obj = seurat_obj,
		func_name = "RunAPOTC",
		id = args$object_id
	)@params

	args$apotc_obj@clone_scale_factor <- args$new_clone_scale_factor

	# FIXME not the best approach - floating point comparison changes shapes
	# TODO use a linear transformation instead: is it s(v_i - c) + c ?
	circlepackClones(
		apotc_obj = args$apotc_obj,
		ORDER = past_params$order_clones,
		try_place = past_params$try_place,
		verbose = args$verbose
	)
}

change_rad_scale <- function(apotc_obj, new_factor) {

	old_factor <- get_rad_scale_factor(apotc_obj)
	conversion_num <- get_clone_scale_factor(apotc_obj) *
		(new_factor - old_factor)

	for (i in seq_len(get_num_clusters(apotc_obj))) {
		curr <- apotc_obj@clusters[[i]]
    	if (isnt_empty(curr)) {
    		apotc_obj@clusters[[i]]$rad <- curr$rad + conversion_num
    	}
	}

	apotc_obj@rad_scale_factor <- new_factor
	apotc_obj
}

recolor_clusters <- function(apotc_obj, recolor_cluster, new_color) {
	for (i in seq_along(recolor_cluster)) {
		apotc_obj@cluster_colors[recolor_cluster[i]] <- new_color[[i]]
	}
	apotc_obj
}

relocate_clusters <- function(apotc_obj, relocate_cluster, relocation_coord) {

	if (is_numeric_pair(relocation_coord)) {
		relocation_coord <- init_list(length(relocate_cluster), relocation_coord)
	}

	new_clusterlists <- get_clusterlists(apotc_obj)

	for (i in seq_along(relocate_cluster)) {
		cl_ind <- relocate_cluster[i]
		if (!is_valid_nonempty_cluster(apotc_obj, cl_ind)) next
		new_clusterlists[[cl_ind]] <- move_cluster(
	    	cluster = new_clusterlists[[cl_ind]],
	    	new_coord = relocation_coord[[i]]
	    )
	}

	setModifiedClusterlists(apotc_obj, new_clusterlists)
}

nudge_clusters <- function(apotc_obj, nudge_cluster, nudge_vector) {

	if (is_numeric_pair(nudge_vector)) {
		nudge_vector <- init_list(length(nudge_cluster), nudge_vector)
	}

	relocate_clusters(
		apotc_obj,
		relocate_cluster = nudge_cluster,
		relocation_coord = operate_on_same_length_lists(
			func = add,
			l1 = nudge_vector,
			l2 = get_centroids(apotc_obj)[nudge_cluster]
		)
	)
}

# label modification stuff - maybe should be in seperate function(s)?

rename_labels <- function(apotc_obj, rename_label, new_label) {
	rename_label <- process_label_input(apotc_obj, rename_label)
	for (i in seq_along(rename_label)) {
		if (!is_valid_nonempty_cluster(apotc_obj, rename_label[i])) next
		apotc_obj@labels[rename_label[i]] <- new_label[i]
	}
	apotc_obj
}

process_label_input <- function(apotc_obj, labels) {
	if (is.numeric(labels)) return(labels)
	if (is.character(labels)) {
		all_labels <- get_labels(apotc_obj)
		labels <- sapply(
			labels, function(x) which(x == all_labels)
		)
	}
	labels
}

relocate_labels <- function(apotc_obj, relocate_label, label_relocation_coord) {
	operate_on_label_locations(
		apotc_obj, relocate_label, label_relocation_coord, function(a, b) b
	)
}

nudge_labels <- function(apotc_obj, nudge_label, label_nudge_vector) {
	operate_on_label_locations(
		apotc_obj, nudge_label, label_nudge_vector, add
	)
}

operate_on_label_locations <- function(apotc_obj, labels, two_d_vectors, func) {

	labels <- process_label_input(apotc_obj, labels)

	if (is_numeric_pair(two_d_vectors)) {
		two_d_vectors <- init_list(length(labels), two_d_vectors)
	}

	for (i in seq_along(labels)) {
		index <- labels[i]
		if (!is_valid_nonempty_cluster(apotc_obj, index)) next
		apotc_obj@label_coords[[index]] <- func(
			apotc_obj@label_coords[[index]],
			two_d_vectors[[i]]
		)
	}

	apotc_obj
}
