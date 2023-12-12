#' @title
#' Directly visualize clonal expansion of a contig combined seurat object
#' 
#' @description
#' `r lifecycle::badge("experimental")`
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
vizAPOTC <- function(
    seurat_obj,
    reduction_base = "umap",
    clonecall = "strict",
    ...,
    extra_filter = NULL,

    clone_scale_factor = "auto",
    rad_scale_factor = 0.95,
    order_clones = TRUE,
    scramble_clones = FALSE,
    try_place = FALSE,
    repulse = TRUE,
    repulsion_threshold = 1,
    repulsion_strength = 1,
    max_repulsion_iter = 20L,

    res = 360,
    linetype = "blank",
    use_default_theme = TRUE,
    retain_axis_scales = FALSE,

    show_labels = FALSE,
    label_size = 5,

    add_size_legend = TRUE,
    legend_sizes = "auto",
    legend_position = "top_left", # can now also be simply a coord
    legend_buffer = 1.5,
    legend_color = "#808080",
    legend_spacing = "auto",
    legend_label = "Clone sizes",
    legend_text_size = 5,

    verbose = TRUE
) {
    # TODO
}