gen_labels <- function(num_clusters) {
	label_vec <- character(num_clusters)
	for (i in seq_len(num_clusters)) {
		label_vec[i] <- paste("C", i, sep = "")
	}
	label_vec
}

# show the labels on the plot
insert_labels <- function(plt, apotc_obj, size) {
	for (i in 1:apotc_obj@num_clusters) {
		if (!isnt_empty(apotc_obj@clusters[[i]])) {
			next
		}

		plt <- plt + ggplot2::annotate(
			"text",
			x = apotc_obj@label_coords[[i]][1],
			y = apotc_obj@label_coords[[i]][2],
			label = apotc_obj@labels[i],
			size = size
		)
	}
	plt
}

# TODO: everything below needs to be fixed with obj_id

modify_names <- function(seurat_obj, apotc_obj, modify_label, to_new_label) {
	label_map <- hash::hash(modify_label, to_new_label)
	for (i in 1:apotc_obj@num_clusters) {
		val <- label_map[[apotc_obj@labels[i]]]
		if (!is.null(val)) {
			apotc_obj@labels[i] <- val
		}
	}
	seurat_obj
}

nudge_label_coords <- function(seurat_obj, modify_label, to_nudge_by) {
	label_map <- hash::hash(modify_label, to_nudge_by)
	for (i in 1:seurat_obj@reduction[['apotc']]@num_clusters) {
		val <- label_map[[seurat_obj@reduction[['apotc']]@labels[i]]]
		if (!is.null(val)) {
			seurat_obj@reduction[['apotc']]@label_coords[i] <- val +
				seurat_obj@reduction[['apotc']]@label_coords[i]
		}
	}
	seurat_obj
}

move_label_coords <- function(seurat_obj, modify_label, to_new_coord) {
	label_map <- hash::hash(modify_label, to_new_coord)
	for (i in 1:seurat_obj@reduction[['apotc']]@num_clusters) {
		val <- label_map[[seurat_obj@reduction[['apotc']]@labels[i]]]
		if (!is.null(val)) {
			seurat_obj@reduction[['apotc']]@label_coords[i] <- val
		}
	}
	seurat_obj
}

#' @title Modify APackOfTheClones cluster labels
#'
#' @description Allows simpler user customization of APackOfTheClones cluster
#' label names and coordinates, and modifies the changes into the input seurat
#' object. Changes will be visible after running `APOTCPlot` on the input
#' seurat object. The length of the lists/vectors must all be the same, as the
#' indicies dictate which labels are changed.
#'
#' @noRd
#' @keywords internal
#'
modifyLabels <- function(
	seurat_obj, modify_label,
	to_new_label = NULL, to_nudge_by = NULL, to_new_coord = NULL
) {
	# needs errors/warnings
	if (!is.null(to_new_label)) {
		seurat_obj <- modify_names(seurat_obj, modify_label, to_new_label)
	}
	if (!is.null(to_nudge_by)) {
		seurat_obj <- nudge_label_coords(seurat_obj, modify_label, to_nudge_by)
	}
	if (!is.null(to_new_coord)) {
		seurat_obj <- move_label_coords(seurat_obj, modify_label, to_new_coord)
	}
	seurat_obj
}
