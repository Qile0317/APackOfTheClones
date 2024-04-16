#' @title Various variations of visualizations of clonal expansion post-RunAPOTC
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Given a seurat object with an 'apotc' (APackOfTheClones) object
#' from running [RunAPOTC], this function will read the information and return
#' a customizable ggplot2 object of the clonal expansion with a circle size
#' legend. If the user is unhappy about certain aspects of the plot, many
#' parameters can be adjusted with the [AdjustAPOTC] function.
#'
#' The specific APackOfTheClones run to be plotted can be identified in two
#' ways: either by inputting the `run_id` associated with the run that was
#' either defined / auto-generated during [RunAPOTC], or by inputting the
#' `reduction_base`, `clonecall`, `extra_filter` and any other keyword arguments
#' that corresponded to the run. Its heavily recommended to use the `run_id`.
#' If none of these parameters are inputted, the function defaults to returning
#' the plot of the latest run.
#'
#' @inheritParams RunAPOTC
#'
#' @param seurat_obj A seurat object that has been integrated with clonotype
#' data and has had a valid run of [RunAPOTC].
#' @param show_shared The output of [getSharedClones] can be inputted here,
#' and the resulting plot will overlay lines between clone circles if that
#' clonotype is common between clusters. Note that the input ***must*** be
#' generated from data in the correct `APackOfTheClones` run, and the behavior
#' is undefined otherwise and will likely error. The next 4 arguments allow for
#' aesthetic customization of these line links.
#' @param only_link Optional integer indicating to only display clone
#' links originating from this cluster if showing shared clones.
#' @param clone_link_width numeric. The width of the lines that connect shared
#' clones. Defaults to `"auto"` which will estimate a reasonable value depending
#' on circle sizes.
#' @param clone_link_color character. The color of the lines that connect shared
#' clones. Defaults to `"blend"` which will use the average colors of the two
#' connected clones. Else, any hex color or valid color string input will work,
#' and the corresponding color will be applied on all links.
#' @param clone_link_alpha numeric. The alpha of the lines that connect shared
#' clones.
#' @param res The number of points on the generated path per full circle. From
#' plot viewers, if circles seem slightly too pixelated, it is recommended to
#' first try to export the plot as an `.svg` before increasing `res` due to
#' increased plotting times from [ggforce::geom_circle].
#' @param linetype The type of outline each circle should have. defaults to
#' `"blank` meaning no outline. More information is in the function
#' documentation of `ggforce::geom_circle`.
#' @param use_default_theme logical that defaults to `TRUE`. If `TRUE`,
#' the resulting plot will have the same theme as the seurat reference reduction
#' plot. Else, the plot will simply have a blank background.
#' @param retain_axis_scales If `TRUE`, approximately maintains the axis scales
#' of the original reduction plot. However, it will only attempt to extend the
#' axes and never shorten. Users are recommended to set this to `TRUE`
#' especially if working with subsetted versions of the clonal data to better
#' preserve the geometric relation to the original dimensional reduction.
#' @param alpha numeric. The alpha of the circles in (0, 1]. Defaults to 1.
#' @param show_labels If `TRUE`, will label each circle cluster at the centroid,
#' defaulting to "C0, C1, ...".
#' @param label_size The text size of labels if shown. Defaults to 5.
#' @param add_size_legend If `TRUE`, adds a legend to the plot visualizing the
#' relative sizes of clones. Note that it is simply an overlay and not a real
#' ggplot2 legend.
#' @param legend_sizes numeric vector. Indicates the circle sizes to be
#' displayed on the legend, and will always be sorted from smallest to greatest.
#' Defaults to `"auto"` which estimate a reasonable range of sizes to display.
#' @param legend_position character or numeric. Can be set to either
#' `"top_left"`, `"top_right"`, `"bottom_left"`, `"bottom_right"` and places the
#' legend roughly in the corresponding position. Otherwise, can be a numeric
#' vector of length 2 indicating the x and y position of the *topmost (smallest)
#' circle* of the legend.
#' @param legend_buffer numeric. Indicates how much to "push" the legend towards
#' the center of the plot from the selected corner. If negative, will push away
#' @param legend_color character. Indicates the hex color of the circles
#' displayed on the legend. Defaults to the hex code for a gray tone
#' @param legend_spacing numeric. Indicates the horizontal distance between each
#' stacked circle on the size legend. Defaults to `"auto"` which will use an
#' estimated value depending on plot size
#' @param legend_label character. The title of the legend, which defaults to
#' `"clone sizes`.
#' @param legend_text_size numeric. The text size of the letters and numbers on
#' the legend
#' @param add_legend_background logical. If `TRUE`, will add a border around the
#' legend and fill the background to be white, overlaying anything else.
#' @param add_legend_centerspace numeric. An additional amount of distance
#' changed between the circle sizes on the left side of the legend and the
#' numbers on the right. Useful to set to around 0.5 (or more / less) when there
#' are particularly large clone sizes that may cover the numbers.
#' @param detail logical. If `FALSE`, will only plot entire clusters as one
#' large circle, which may be useful in cases where there are a high number
#' of clones resulting in a large number of circles on the resulting ggplot,
#' which has increased plotting times, and certain aspects of the plot needs
#' to be finely adjusted with [AdjustAPOTC] or simply inspected. This should
#' not be set to `FALSE` for the actual clonal expansion plot.
#'
#' @return A ggplot object of the APackOfTheClones clonal expansion plot of the
#' seurat object. There is an additional 10th element in the object named
#' `"APackOfTheClones"` used by other functions in this package and shouldn't
#' interfere with any other ggplot functionality. (As far as currently known)
#' @export
#'
#' @seealso [AdjustAPOTC]
#'
#' @examples
#' data("combined_pbmc")
#'
#' combined_pbmc <- RunAPOTC(
#'     combined_pbmc, run_id = "run1", verbose = FALSE
#' )
#'
#' # plotting with default arguments will plot the latest "run1"
#' clonal_packing_plot <- APOTCPlot(combined_pbmc)
#'
APOTCPlot <- function(
	seurat_obj,
	reduction_base = NULL,
	clonecall = NULL,
	...,
	extra_filter = NULL,
	run_id = NULL,

	show_shared = NULL,
	only_link = NULL,
	clone_link_width = "auto",
	clone_link_color = "black",
	clone_link_alpha = 0.5,

	res = 360L,
	linetype = "blank",
	use_default_theme = TRUE,
	retain_axis_scales = FALSE,
	alpha = 1,

	show_labels = FALSE,
	label_size = 5,

	add_size_legend = TRUE,
	legend_sizes = "auto",
	legend_position = "auto",
	legend_buffer = 0.2,
	legend_color = "#808080",
	legend_spacing = "auto",
	legend_label = "Clone sizes",
	legend_text_size = 5,
	add_legend_background = TRUE,
	add_legend_centerspace = 0,

	detail = TRUE,

	verbose = TRUE
) {
	# handle varargs, run_id, and typecheck
	varargs_list <- list(...)
	args <- environment()
	args$run_id <- infer_object_id_if_needed(args, varargs_list)
	APOTCPlot_error_handler(args)

	# get the apotc object
	apotc_obj <- getApotcData(seurat_obj, args$run_id)

	# initialize plot
	result_plot <- create_initial_apotc_plot(
		apotc_obj, res, linetype, alpha, detail
	)
	result_plot_dimensions <- get_apotc_plot_dims(apotc_obj)

	#set theme
	if (use_default_theme) {
		result_plot <- add_default_theme(
			plt = result_plot,
			reduction = get_reduction_base(apotc_obj)
		)
	} else {
		result_plot <- result_plot + ggplot2::theme_void()
	}

	# retain axis scales on the resulting plot.
	if (retain_axis_scales) {
		result_plot_dimensions <- get_retain_scale_dims(
			seurat_obj,
			reduction = get_reduction_base(apotc_obj),
			ball_pack_plt = result_plot,
			plot_dims = result_plot_dimensions
		)

		result_plot <- result_plot + ggplot2::expand_limits(
			x = get_xr(result_plot_dimensions),
			y = get_yr(result_plot_dimensions)
		)
	}

	if (isnt_empty(show_shared)) {

		# check only_link indexing
		if (!is.null(only_link) &&
			!is_valid_nonempty_cluster(apotc_obj, only_link)) {
			warning(call. = FALSE,
				"* The cluster at index `only_link` = ", only_link,
				" is empty or isn't between 1 ~ ", get_num_clusters(apotc_obj)
			)
			only_link <- NULL
		}

		result_plot <- overlay_shared_clone_links(
			apotc_obj = apotc_obj,
			shared_clones = show_shared,
			result_plot = result_plot,
			only_cluster = only_link,
			link_color_mode = clone_link_color,
			link_width = clone_link_width,
			link_alpha = clone_link_alpha,
			verbose = verbose
			# TODO other params in the future
		)
	}

	if (show_labels) {
		result_plot <- insert_labels(result_plot, apotc_obj, label_size)
	}

	if (add_size_legend) {
		result_plot <- insert_legend(
			plt = result_plot,
			plt_dims = result_plot_dimensions,
			apotc_obj = apotc_obj,
			sizes = legend_sizes,
			pos = legend_position,
			buffer = legend_buffer,
			additional_middle_spacing = add_legend_centerspace,
			color = legend_color,
			n = res,
			spacing = legend_spacing,
			legend_label = legend_label,
			legend_textsize = legend_text_size,
			do_add_legend_border = add_legend_background,
			linetype = linetype
		)
	}

	result_plot <- ApotcGGPlot(result_plot, apotc_obj)

	if (verbose) message("* generated ggplot object")
	result_plot
}

APOTCPlot_error_handler <- function(args) {

	check_apotc_identifiers(args)
	check_filtering_conditions(args)

	# check object_id validity
	if (!containsApotcRun(args$seurat_obj, args$run_id)) {
		stop(call. = FALSE, paste(
			"APackOfTheClones object with id", args$run_id,
			"does not exist in the seurat object"
		))
	}

	# typecheck clone link args
	typecheck(args$show_shared, is_output_of_getSharedClones, is.null)
	typecheck(args$only_link, is_an_integer, is.null)
	if (!should_estimate(args$clone_link_width))
		typecheck(args$clone_link_width, is_a_positive_numeric)
	typecheck(args$clone_link_color, is_a_character)
	typecheck(args$clone_link_alpha, is_a_numeric)

	# typecheck visualization args
	typecheck(args$res, is_an_integer)
	typecheck(args$linetype, is_a_character)
	typecheck(args$use_default_theme, is_a_logical)
	typecheck(args$retain_axis_scales, is_a_logical)
	typecheck(args$show_labels, is_a_logical)
	typecheck(args$label_size, is_a_positive_numeric)

	# check legend args
	typecheck(args$add_size_legend, is_a_logical)
	check_legend_params(args)

	typecheck(args$detail, is_a_logical)

}

# helpers for getting plot dimensions quickly

get_apotc_plot_dims <- function(apotc_obj) {
	apotc_obj %>%
		get_plottable_df_with_color() %>%
		get_apotc_plot_dims_from_df()
}

get_apotc_plot_dims_from_df <- function(plot_dataframe) {
	plot_dataframe %>%
		subset_to_only_edge_circles() %>%
		plot_clusters() %>%
		get_plot_dims()
}

subset_to_only_edge_circles <- function(apotc_plot_dataframe) {
	apotc_plot_dataframe[unique(rcppGetEdgeCircleindices(apotc_plot_dataframe)), ]
}

# Produce modified ggplot object of an APackOfTheClones plot with an extra slot
# hack fix - the clone size slot stores the autogenerated legend sizes
ApotcGGPlot <- function(ggplot_obj, apotc_obj) {
	apotc_obj@clusters <- list()
	apotc_obj@clone_sizes <- list(estimate_legend_sizes(apotc_obj))
	ggplot_obj$APackOfTheClones <- apotc_obj
	ggplot_obj
}

isApotcGGPlot <- function(ggplot_obj) {
	inherits(ggplot_obj, "ggplot") && !is.null(ggplot_obj$APackOfTheClones)
}

# getter
get_apotcdata <- function(apotc_ggplot_obj) {
	apotc_ggplot_obj$APackOfTheClones
}

# based on the hack fix - get the estimated legend sizes in a apot ggplot
get_estimated_legend_sizes <- function(apotc_ggplot_obj) {
	get_raw_clone_sizes(get_apotcdata(apotc_ggplot_obj))[[1]]
}

