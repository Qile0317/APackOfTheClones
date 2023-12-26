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
#' @param seurat_obj Object of class `Seurat` which has had `integrate_tcr` and
#' `RunAPOTC` functions ran on.
#' @param modify_label character (vector) of labels to modify (to rename
#' and/or move)
#' @param to_new_label character (vector) of new labels to convert the labels in
#' `modify_label` to. The order MUST be the exact same as in modify_label
#' @param `to_nudge_by` list of (or one single) numeric vector of length two
#' indicating the (x,y) amounts to *shift* the corresponding label by in the plot
#' @param `to_new_coord` list of (or one single) numeric vector of length two
#' indicating the (x,y) amounts to *move* the corresponding label to in the plot
#'
#' @return A modified version of the input `seurat_obj` where
#' `@reduction[['apotc']]@labels` and `@reduction[['apotc']]@label_coords` were
#' accordingly modified
#'
#' @export
#'
#' @examples
#' # TODO all of these need to be redone
#' # rename "C0" and move it to x = 1, y = 2
#' seurat_obj <- modifyLabels(
#'     seurat_obj, modify_label = "C0",
#'     to_new_label = "CD4+ CTL", to_new_coord = c(1, 2)
#' )
#'
#' # rename multiple labels
#' seurat_obj <- modifyLabels(
#'     seurat_obj, modify_label = c("C4", "C9", "C2"),
#'     to_new_label = c("T cm", "T cm 2", "T fh")
#' )
#'
#' # shift multiple labels
#' seurat_obj <- modifyLabels(
#'     seurat_obj, modify_label = c("T cm 2", "T fh", "C12"),
#'     to_nudge_by = list(c(1,1), c(4,3), c(0.2,-2))
#' )
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
