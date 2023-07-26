gen_labels <- function(num_clusters) {
	label_vec <- vector('character', num_clusters)
	for (i in 0:num_clusters-1) {
		label_vec[i] <- paste("C", i, sep = "")
	}
	label_vec
}

insert_labels <- function(plt, seurat_obj, size) {
	for (i in 1:seurat_obj@reductions[['apotc']]@num_clusters) {
		if (!isnt_empty(seurat_obj@reductions[['apotc']]@clusters[[i]])) {
			next
		}

		plt <- plt + ggplot2::annotate(
			"text",
			x = seurat_obj@reductions[['apotc']]@label_coords[1],
			y = seurat_obj@reductions[['apotc']]@label_coords[2],
			label = seurat_obj@reductions[['apotc']]@labels[i],
			size = size
		)
	}
	plt
}

ModifyLabels <- function(seurat_obj) {

}
