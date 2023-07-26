add_default_theme <- function(plt, reduction) {
	label_hashmap <- hash::hash(
		c("umap", "tsne", "pca"), c("UMAP", "tSNE", "PC")
	)
	label <- label_hashmap[[reduction]]

	plt +
		ggplot2::theme_classic() +
		ggplot2::xlab(paste(label, 1, sep = "_")) +
		ggplot2::ylab(paste(label, 2, sep = "_")) +
		ggplot2::ggtitle("Sizes of clones within each cluster")
}

#' Create clonal expansion plot after RunAPOTC()
#'
#' Given a seurat object with an 'apotc' (APackOfTheClones) object from running
#' RunAPOTC(), this function will read the information and return a customizable
#' ggplot of the clonal expansion with a circle size legend.
#'
#'
APOTCPlot <- function(
	seurat_obj,
	res = 360,
	linetype = "blank",
	use_default_theme = TRUE,
	show_origin = FALSE,
	retain_axis_scales = FALSE,

	show_labels = FALSE,
	label_size = 5,

	add_size_legend = TRUE,
	legend_sizes = c(1, 5, 50),
	legend_position = "top_left", # can now also be simply a coord
	legend_buffer = 1.5,
	legend_color = "#808080",
	legend_spacing = "auto",
	legend_label = "Clone sizes",
	legend_text_size = 5
) {
	if (is.null(seurat_obj@reductions[["apotc"]])) {
		stop("Please do an APackOfTheClones run on the seurat object first with RunAPOTC")
	}

	# convert clusterlists to dataframe and add colors
	clusterlists <- seurat_obj@reductions[["apotc"]]@clusters
	clusterlists <- df_full_join(clusterlists)
	clusterlists <- extract_and_add_colors(seurat_obj, clusterlists) # will change

	result_plot <- plot_clusters(
		clusterlists, n = res, linetype = linetype,
		title = "Sizes of Clones Within Each Cluster", haslegend = FALSE,
		void = FALSE, origin = show_origin
	)

	#set theme
	if (use_default_theme) {
		result_plot <- add_default_theme(
			result_plot, seurat_obj@reductions[["apotc"]]@reduction_base
		)
	} else {
		result_plot <- result_plot + ggplot2::theme_void()
	}

	# retain axis scales on the resulting plot. The function sucks tho
	if (retain_axis_scales) {
		result_plot <- suppressMessages(invisible(
			retain_scale(seurat_obj, result_plot)
		))
	}

	if (add_size_legend) {
		rad_scale <- seurat_obj@reductions[["apotc"]]@rad_scale_factor
		clone_scale <- seurat_obj@reductions[["apotc"]]@clone_scale_factor

		result_plot <- insert_legend(
			plt = result_plot, circ_scale_factor = clone_scale,
			rad_decrease = convert_to_rad_decrease(rad_scale, clone_scale),
			sizes = legend_sizes, pos = legend_position, buffer = legend_buffer,
			color = legend_color, n = res, spacing = legend_spacing,
			legend_label = legend_label, legend_textsize = legend_text_size
		)
	}

	if (show_labels) {
		result_plot <- insert_labels(result_plot, seurat_obj, label_size)
	}

	result_plot
}
