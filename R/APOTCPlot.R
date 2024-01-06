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
#'
#' @return A ggplot object of the APackOfTheClones clonal expansion plot of the
#' seurat object
#'
#' @seealso [AdjustAPOTC]
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' suppressPackageStartupMessages(library(APackOfTheClones))
#'
#' # 2 arbitrary APackOfTheClones runs named "run1" and "run2"
#' combined_pbmc <- RunAPOTC(get(data("combined_pbmc")), run_id = "run1")
#' combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "run2", try_place = TRUE)
#'
#' # plotting with default arguments will plot the "run2" results as its latest
#' APOTCPlot(combined_pbmc)
#'
#' # otherwise, one should plot with the run_id argument for other runs
#' APOTCPlot(combined_pbmc, run_id = "run2")
#'
APOTCPlot <- function(
	seurat_obj,
	reduction_base = NULL,
	clonecall = NULL,
	...,
	extra_filter = NULL,
	run_id = NULL,

	res = 360L,
	linetype = "blank",
	use_default_theme = TRUE,
	retain_axis_scales = FALSE,
	#alpha = 1,

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
	add_legend_background = TRUE
) {

	args <- hash::hash(as.list(environment()))
	APOTCPlot_error_handler(args)

	apotc_obj <- getApotcData(
		seurat_obj, infer_object_id_if_needed(args, varargs_list = list(...))
	)

	result_plot <- plot_clusters(
		clusters = get_plottable_df_with_color(apotc_obj),
		n = res,
		linetype = linetype#,
		#alpha=alpha
	)

	#set theme
	if (use_default_theme) {
		result_plot <- add_default_theme(
			plt = result_plot,
			reduction = get_reduction_base(apotc_obj)
		)
	} else {
		result_plot <- result_plot + ggplot2::theme_void()
	}

	# retain axis scales on the resulting plot. The function sucks tho
	if (retain_axis_scales) {
		result_plot <- suppressMessages(invisible(
			retain_scale(seurat_obj, get_reduction_base(apotc_obj), result_plot)
		))
	}

	# TODO clonal link computation here

	if (show_labels) {
		result_plot <- insert_labels(result_plot, apotc_obj, label_size)
	}

	if (add_size_legend) {
		result_plot <- insert_legend(
			plt = result_plot,
			plt_dims = get_plot_dims(result_plot),
			apotc_obj = apotc_obj,
			sizes = legend_sizes,
			pos = legend_position,
			buffer = legend_buffer,
			color = legend_color,
			n = res,
			spacing = legend_spacing,
			legend_label = legend_label,
			legend_textsize = legend_text_size,
			do_add_legend_border = add_legend_background
		)
	}
	
	result_plot
}

APOTCPlot_error_handler <- function(args) {
	if (!is_seurat_object(args[["seurat_obj"]])) {
        stop("`seurat_obj` must be a Seurat object.")
    }

	if (!is.null(args$reduction_base) && !is_a_character(args$reduction_base)) {
		stop("`reduction_base` must be a character of length 1.")
	}
    
	if (!is.null(args$clonecall) && !is_a_character(args$clonecall)) {
		stop("`clonecall` must be a character of length 1.")
	}

	if (!is.null(args$extra_filter) && !is_a_character(args$extra_filter)) {
		stop("`extra_filter` must be a character or NULL of length 1.")
	}

    if (!is_an_integer(args$res)) {
        stop("`res` must be an integer value of length 1.")
    }

    if (!is_a_character(args$linetype)) {
        stop("`linetype` must be a character of length 1.")
    }

    if (!is_a_logical(args$use_default_theme)) {
        stop("`use_default_theme` must be a logical value of length 1.")
    }

    if (!is_a_logical(args$retain_axis_scales)) {
        stop("`retain_axis_scales` must be a logical value of length 1.")
    }

    if (!is_a_logical(args$show_labels)) {
        stop("`show_labels` must be a logical value of length 1.")
    }

    if (!is_a_numeric(args$label_size)) {
        stop("`label_size` must be a numeric value of length 1.")
    }

    if (!is_a_logical(args$add_size_legend)) {
        stop("`add_size_legend` must be a logical value of length 1.")
    }

	# check object_id validity
	if (!containsApotcRun(args$seurat_obj, args$run_id)) {
		stop(paste(
			"APackOfTheClones object with id", args$run_id,
			"does not exist in the seurat object"
		))
	}

    # TODO: Add more specific checks for other parameters
}
