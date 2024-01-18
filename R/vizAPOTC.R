#' @title
#' Directly visualize clonal expansion of a combined seurat object
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function combines the functionality of both `RunAPOTC` and `APOTCPlot`.
#' Given a Seurat object, it first runs the APackOfTheClones method ([RunAPOTC])
#' to compute clonal expansion information, and then generates a customizable
#' ggplot2 object of the clonal expansion plot with a circle size legend
#' ([APOTCPlot]).
#'
#' @inheritParams RunAPOTC
#' @inheritParams APOTCPlot
#'
#' @inherit APOTCPlot return
#' @export
#'
#' @examples
#' data("combined_pbmc")
#'
#' # plot with default parameters
#' vizAPOTC(combined_pbmc, verbose = FALSE)
#'
#' # use arguments from RunAPOTC and APOTCPlot
#' vizAPOTC(
#'     combined_pbmc, try_place = TRUE, show_labels = TRUE, verbose = FALSE
#' )
#'
vizAPOTC <- function(
    seurat_obj,
    reduction_base = "umap",
    clonecall = "strict",
    ...,
    extra_filter = NULL,

    clone_scale_factor = "auto",
    rad_scale_factor = 0.95,
    order_clones = TRUE,
    try_place = FALSE,
    
    repulse = TRUE,
    repulsion_threshold = 1,
    repulsion_strength = 1,
    max_repulsion_iter = 20L,

    show_shared_clones = FALSE,

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
    add_legend_background = TRUE,

    verbose = TRUE
) {
    seurat_obj <- RunAPOTC(
        seurat_obj,
        reduction_base = reduction_base,
        clonecall = clonecall,
        ...,
        extra_filter = extra_filter,
        run_id = "vizAPOTC",

        clone_scale_factor = clone_scale_factor,
        rad_scale_factor = rad_scale_factor,
        order_clones = order_clones,
        try_place = try_place,
        repulse = repulse,
        repulsion_threshold = repulsion_threshold,
        repulsion_strength = repulsion_strength,
        max_repulsion_iter = max_repulsion_iter,
        override = TRUE,
        verbose = verbose
    )

    if (verbose) message("Plotting...\n")

    APOTCPlot(
        seurat_obj,
        run_id = "vizAPOTC",

        show_shared_clones = show_shared_clones,

        res = res,
        linetype = linetype,
        use_default_theme = use_default_theme,
        retain_axis_scales = retain_axis_scales,
        
        show_labels = show_labels,
        label_size = label_size,

        add_size_legend = add_size_legend,
        legend_sizes = legend_sizes,
        legend_position = legend_position,
        legend_buffer = legend_buffer,
        legend_color = legend_color,
        legend_spacing = legend_spacing,
        legend_label = legend_label,
        legend_text_size = legend_text_size,
        add_legend_background = add_legend_background
    )
}
