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

APOTCPlot <- function(
	seurat_obj,
	res = 360,
	linetype = "blank",
	use_default_theme = TRUE,
	show_origin = FALSE,
	retain_axis_scales = FALSE,
	add_size_legend = TRUE,
	legend_sizes = c(1, 5, 50),
	legend_position = "top_left", # can now also be simply a coord
	legend_buffer = 1.5,
	legend_color = "#808080",
	legend_spacing = 0.4
) {
	if (is.null(seurat_obj@reductions[["apotc"]])) {
		stop("Please do an APackOfTheClones run on the seurat object first")
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
		return(insert_legend(
			plt = result_plot,
			circ_scale_factor = seurat_obj@reductions[["apotc"]]@clone_scale_factor,
			rad_scale_factor = seurat_obj@reductions[["apotc"]]@rad_scale_factor,
			sizes = legend_sizes, pos = legend_position, buffer = legend_buffer,
			color = legend_color, n = res, spacing = legend_spacing
		))
	}
	result_plot
}
