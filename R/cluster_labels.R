gen_labels <- function(num_clusters) {
	label_vec <- character(num_clusters)
	for (i in seq_len(num_clusters)) {
		label_vec[i] <- paste("C", i, sep = "")
	}
	label_vec
}

# show the labels on the plot
insert_labels <- function(plt, apotc_obj, size) {
	for (i in seq_len(get_num_clusters(apotc_obj))) {
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
