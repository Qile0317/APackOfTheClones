# the following code is copied and modified from ggplot2 under the MIT license

apotc_rename_aes <- function(x) {
    names(x) <- ggplot2::standardise_aes_names(names(x))
    duplicated_names <- names(x)[duplicated(names(x))]
    if (length(duplicated_names) > 0L) {
        return(NULL)
    }
    x
}

apotc_aes_string <- function(x, y, ...) {
    mapping <- list(...)
    if (!missing(x)) mapping["x"] <- list(x)
    if (!missing(y)) mapping["y"] <- list(y)

    mapping <- lapply(mapping, function(x) {
        if (is.character(x)) {
            x <- rlang::parse_expr(x)
        }
    })
    structure(apotc_rename_aes(mapping), class = "uneval")
}

#' @title Get the xmin, xmax, ymin, ymax of a ggplot object
#' @return list(xr = c(xmin, xmax), yr = c(ymin, ymax))
#' @noRd
get_plot_dims <- function(plt) {
    built_plt_layout <- ggplot2::ggplot_build(plt)$layout
    list(
        xr = built_plt_layout$panel_scales_x[[1]]$range$range,
        yr = built_plt_layout$panel_scales_y[[1]]$range$range
    )
}

get_xr <- function(p) {
    if (ggplot2::is_ggplot(p)) {
        return(ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$range$range)
    }
    p[[1]]
}

get_yr <- function(p) {
    if (ggplot2::is_ggplot(p)) {
        return(ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    }
    p[[2]]
}

name_latest_layer <- function(plt, new_name) {
    if (is.null(names(plt$layers))) {
        names(plt$layers) <- rep("", length(plt$layers))
    }
    names(plt$layers)[length(plt$layers)] <- new_name
    plt
}

remove_ggplot_layers <- function(ggplot_obj, layer_indices) {
    ggplot_obj$layers[layer_indices] <- NULL
    ggplot_obj
}

get_ggplot_data <- function(x) x$data

set_ggplot_data <- function(ggplot_obj, new_data) {
    ggplot_obj$data <- new_data
    ggplot_obj
}
