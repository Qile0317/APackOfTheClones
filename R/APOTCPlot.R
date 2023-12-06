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

#' @title Create clonal expansion plot after RunAPOTC()
#'
#' @description Given a seurat object with an 'apotc' (APackOfTheClones) object
#' from running `RunAPOTC()`, this function will read the information and return
#' a customizable ggplot object of the clonal expansion with a circle size legend.
#' If the user is unhappy about certain aspects of the plot, some parameters can
#' be adjusted with the `AdjustAPOTC` function.
#'
#' @param seurat_obj A seurat object that has been integrated with clonotype
#' data and has had `RunAPOTC` ran on it.
#' @param res The number of points on the generated path per full circle. From
#' plot viewers, if circles seem slightly too pixelated, it is recommended to
#' first try to export the plot as an `.svg` before increasing `res` due to
#' increased plotting times from `ggforce::geom_circle`.
#' @param linetype The type of outline each circle should have. defaults to
#' `"blank` meaning no outline. More information is in the function
#' documentation of `ggforce::geom_circle`.
#' @param use_default_theme logical that defaults to `TRUE`. If `TRUE`,
#' the resulting plot will have the same theme as the seurat reference reduction
#' plot. Else, the plot will simply have a blank background.
#' @param retain_axis_scales If `TRUE`, approximately maintains the axis scales
#' of the original UMAP. However, it will only attempt to extend the axes and
#' never shorten.
#' @param show_labels If `TRUE`, will label each circle cluster at the centroid,
#' defaulting to "C0, C1, ...". labels and their coordinates can also be
#' manually edited with the `ModifyLabels` function.
#' @param label_size The text size of labels if shown. Defaults to 5.
#' @param add_size_legend If `TRUE`, adds a legend to the plot visualizing the
#' relative sizes of clones. Note that it is simply an overlay and not a real
#' ggplot2 legend.
#' @param legend_sizes numeric vector. Indicates the circle sizes to be
#' displayed on the legend and defaults to `c(1, 5, 10)`.
#' @param legend_position character or numeric. Can be set to either
#' `"top_left"`, `"top_right"`, `"bottom_left"`, `"bottom_right"` and places the
#' legend roughly in the corresponding position. Otherwise, can be a numeric
#' vector of length 2 indicating the x and y position of the "top-center" of the
#' legend
#' @param legend_buffer numeric. Indicates how much to "push" the legend towards the center of the plot from the selected corner. If negative, will push away
#' @param legend_color character. Indicates the hex color of the circles displayed on the legend. Defaults to the hex code for gray
#' @param legend_spacing numeric. Indicates the horizontal distance between each
#' stacked circle on the size legend. Defaults to `"auto"` which will use an
#' estimated value depending on plot size
#' @param legend_label character. The title of the legend, which defaults to
#' `"clone sizes`.
#' @param legend_text_size numeric. The text size of the letters and numbers on
#' the legend
#'
#' @return A ggplot object of the APackOfTheClones clonal expansion plot of the
#' seurat object
#'
#' @seealso [AdjustAPOTC()]
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' suppressPackageStartupMessages(library(APackOfTheClones))
#' data("mini_clonotype_data","mini_seurat_obj")
#'
#' # first the APOTC pipeline has to be run
#' pbmc <- integrate_tcr(mini_seurat_obj, mini_clonotype_data, verbose = FALSE)
#' pbmc <- RunAPOTC(pbmc, verbose = FALSE)
#'
#' # generate the default plot with APOTCPlot
#' APOTCPlot(pbmc)
#'
#' # if plotting of the same object with different customizations
#' APOTCPlot(pbmc, use_default_theme = FALSE, show_labels = TRUE)
#'
APOTCPlot <- function(
	seurat_obj,
	res = 360,
	linetype = "blank",
	use_default_theme = TRUE,
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
		void = FALSE, origin = FALSE
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
