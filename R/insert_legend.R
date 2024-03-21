# script to manage the circle size legend overlay

#' @title Remove current APackOfTheClones legend
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Removes the clone size legend on an APackOfTheClones plot, if
#' one is present. Will preserve any additional ggplot layers.
#'
#' @inheritParams overlayLegend
#' @inheritParams APOTCPlot
#' 
#' @inherit APOTCPlot return
#' @export
#'
#' @seealso [overlayLegend]
#'
#' @examples
#' # create an APackOfTheClones plot with a legend
#' apotc_plot <- vizAPOTC(
#'     get(data("combined_pbmc")),
#'     add_size_legend = TRUE,
#'     verbose = FALSE
#' )
#'
#' # remove the legend
#' apotc_plot <- removeLegend(apotc_plot)
#' apotc_plot
#'
removeLegend <- function(apotc_ggplot) {
    check_is_apotc_ggplot(apotc_ggplot)
    if (!has_legend(apotc_ggplot)) return(apotc_ggplot)
    remove_ggplot_layers(apotc_ggplot, get_legend_layer_indicies(apotc_ggplot))
}

#' @title overlay a clone size legend on an APackOfTheClones plot
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function has most of the parameters related to legend in
#' [APOTCPlot], and can plot a new / override the current legend.
#' However, it is very important that the input plot to the function
#' is a plot generated solely by [APOTCPlot] or [vizAPOTC] due to it
#' being a custom ggplot object. It will not override or erase any
#' additional layers that the user/other functions have added. To just
#' remove the legend, see [removeLegend].
#'
#' @param apotc_ggplot a ggplot object that is the output of [APOTCPlot] or
#' [vizAPOTC]
#' @inheritParams APOTCPlot
#'
#' @details
#' The size legend on APackOfTheClones plots are simply a collection of
#' annotation layers of the `rect`, `text`, `circle` geoms. Therefore it
#' isn't quite a ggplot legend. In the actual ggplot object, the `$layers`
#' element should all be named with an empty character `""`, and those
#' that comprise the annotation layers of the legend should be named whatever
#' `APackOfTheClones:::.ApotcLegendLayerName` is. Note that this is simply
#' an implementation detail that the user should not interfere with.
#'
#' @inherit APOTCPlot return
#' @export
#'
#' @seealso [removeLegend]
#'
#' @examples
#' # create a plot with a legend
#' apotc_plot <- vizAPOTC(get(data("combined_pbmc")), verbose = FALSE)
#'
#' # reposition the legend to top right
#' overlayLegend(apotc_plot, legend_position = "top right")
#'
#' # use different sizes and label - may be nice to use the pipe
#' apotc_plot |> overlayLegend(
#'     legend_sizes = c(1, 3, 7, 9),
#'     legend_label = "odd sizes"
#' )
#'
overlayLegend <- function(
    apotc_ggplot,
    legend_sizes = "auto",
	legend_position = "auto",
	legend_buffer = 0.2,
	legend_color = "#808080",
	legend_spacing = "auto",
	legend_label = "Clone sizes",
	legend_text_size = 5,
	add_legend_background = TRUE,
	add_legend_centerspace = 0,
    linetype = "blank",
    res = 360L
) {
    
    overlayLegend_error_handler()

    layers_after_legend_present <- FALSE

    if (has_legend(apotc_ggplot)) {
        layers_after_legend <- get_layers_after_legend(apotc_ggplot)
        layers_after_legend_present <- !is.null(layers_after_legend)
        apotc_ggplot <- remove_legend_and_layers_after(apotc_ggplot)
    }

    if (should_estimate(legend_sizes)) {
        legend_sizes <- get_estimated_legend_sizes(apotc_ggplot)
    }

    apotc_ggplot <- insert_legend(
        plt = apotc_ggplot,
        plt_dims = get_apotc_plot_dims_from_df(apotc_ggplot$data),
        apotc_obj = get_apotcdata(apotc_ggplot),
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

    if (layers_after_legend_present) {
        apotc_ggplot$layers <- append(apotc_ggplot$layers, layers_after_legend)
    }

    apotc_ggplot
}

overlayLegend_error_handler <- function() {
    args <- get_parent_func_args()
    check_legend_params(args)
    check_is_apotc_ggplot(args$apotc_ggplot)
    typecheck(args$res, is_a_positive_integer)
}

check_legend_params <- function(args) {

    if (!should_estimate(args$legend_sizes))
        typecheck(args$legend_sizes, is_positive_numeric)

    if (!(
        is_a_character(args$legend_position) ||
            is_numeric_pair(args$legend_position)
    )) {
        stop(call. = FALSE,
            "`legend_position` must be a numeric pair or specific characters. ",
            "See function level documentation for details on this argument."
        )
    }

    typecheck(args$legend_buffer, is_a_numeric)
    typecheck(args$legend_color, is_a_character)

    if (!(
        is_a_numeric(args$legend_spacing) ||
            should_estimate(args$legend_spacing)
    )) {
        stop(call. = FALSE, "`legend_spacing` must be a numeric.")
    }

    typecheck(args$legend_label, is_a_character)
    typecheck(args$legend_text_size, is_a_positive_numeric)
    typecheck(args$add_legend_background, is_a_logical)
    typecheck(args$add_legend_centerspace, is_a_numeric)
}

# internal global constant for the name of each legend layer
.ApotcLegendLayerName <- "ApotcLegendLayer"
utils::globalVariables(".ApotcLegendLayerName")

name_latest_legend_layer <- function(plt) {
    plt %>% name_latest_layer(.ApotcLegendLayerName)
}

has_legend <- function(apotc_ggplot) {
    any(names(apotc_ggplot$layers) == .ApotcLegendLayerName)
}

has_layers_after_legend <- function(apotc_ggplot) {
    length(apotc_ggplot$layers) > get_last_legend_layer_index(apotc_ggplot)
}

get_layers_after_legend <- function(apotc_ggplot) {
    if (!has_layers_after_legend(apotc_ggplot)) return(NULL)
    layers <- apotc_ggplot$layers
    layers[(get_last_legend_layer_index(apotc_ggplot) + 1):length(layers)]
}

get_legend_layer_indicies <- function(apotc_ggplot) {
    which(names(apotc_ggplot$layers) == .ApotcLegendLayerName)
}

get_first_legend_layer_index <- function(apotc_ggplot) {
    get_legend_layer_indicies(apotc_ggplot)[1]
}

get_last_legend_layer_index <- function(apotc_ggplot) {
    getlast(get_legend_layer_indicies(apotc_ggplot))
}

remove_legend_and_layers_after <- function(apotc_ggplot) {
    num_layers <- length(apotc_ggplot$layers)
    removal_indicies <- get_first_legend_layer_index(apotc_ggplot):num_layers
    apotc_ggplot$layers[removal_indicies] <- NULL
    apotc_ggplot
}

insert_legend <- function(
    plt,
    plt_dims,
    apotc_obj,
    sizes,
    pos,
    buffer,
    additional_middle_spacing,
    color = "#808080",
    n = 360,
    spacing = "auto",
    legend_label = "Clone sizes",
    legend_textsize = 5,
    do_add_legend_border = FALSE,
    linetype = "blank"
) {

    # setup relevant variables

    rad_decrease <- get_rad_decrease(apotc_obj)
    sizes <- get_processed_legend_sizes(apotc_obj, sizes)
    pos <- correct_legend_coord_if_str(pos)
    spacing <- process_legend_spacing(spacing, plt_dims, rad_decrease)

    # calculate relevant legend plotting data

    unpositioned_legend_df <- gen_unpositioned_legend_df(
        legend_sizes = sizes,
        spacing = spacing,
        circ_scale_factor = apotc_obj@clone_scale_factor,
        rad_decrease = rad_decrease,
        color = color,
        additional_middle_spacing = additional_middle_spacing
    )

    legend_dims <- get_legend_dims(unpositioned_legend_df)

    if (should_estimate(pos)) {
        legend_df <- estimate_best_legend_df(
            plt_dims, plt$data, unpositioned_legend_df, legend_dims, buffer
        )
    } else {
        legend_df <- generate_legend_df(
            pos, unpositioned_legend_df, legend_dims, plt_dims, buffer
        )
    }

    # plotting

    # add the legend label on top
    label_coord <- get_legend_title_coord(legend_df, legend_dims, spacing)
    plt <- (plt + ggplot2::annotate(
        "text", x = label_coord[1], y = label_coord[2],
        label = legend_label, size = legend_textsize
    )) %>% name_latest_legend_layer()

    # add the background
    if (do_add_legend_border) {
        plt <- add_legend_backing(plt, legend_df) %>% name_latest_legend_layer()
    }
    
    # add the side number labels
    plt <- (plt + ggplot2::annotate(
        "text", x = legend_df[, "label_x"], y = legend_df[, "y"],
        label = legend_df[, "labels"], size = legend_textsize
    )) %>% name_latest_legend_layer()

    # add the circles and return
    (plt + ggforce::geom_circle(
        data = legend_df,
        mapping = apotc_aes_string(
            x0 = "circle_x",
            y0 = "y",
            r = "rad",
            fill = "color"
        ),
        linetype = linetype,
        n = n
    )) %>% name_latest_legend_layer()
}

get_processed_legend_sizes <- function(apotc_obj, s) {
    if (is.numeric(s)) return(sort(unique(c(1, s))))
    if (should_estimate(s)) return(estimate_legend_sizes(apotc_obj))
    s
}

estimate_legend_sizes <- function(apotc_obj) {
    sizes <- unlist(get_raw_clone_sizes(apotc_obj))
    sort(unique(round(c(
        1,
        stats::median(sizes), mean(sizes), mean(unique(sizes)),
        max(sizes)
    ))))
}

correct_legend_coord_if_str <- function(pos) {
    if (is.numeric(pos)) {
        return(pos)
    }

    pos <- strip_and_lower(pos)
    pos <- switch(
        pos,
        "auto" = "auto",
        "topleft" = "top_left",
        "topright" = "top_right",
        "bottomleft" = "bottom_left",
        "bottomright" = "bottom_right",
        pos
    )

    if (should_estimate(pos)) return(pos)

    user_attempt_correction(
        s = pos,
        strset = c("top_left", "top_right", "bottom_left", "bottom_right"),
        stop_msg_start = "invalid legend coordinate string"
    )
}

process_legend_spacing <- function(spacing, plt_dims, rad_decrease) {
    if (should_estimate(spacing)) {
        return(calculate_legend_spacing(spacing, plt_dims, rad_decrease))
    }
    spacing
}

calculate_legend_spacing <- function(
    spacing, plt_dims, rad_decrease, portion = 0.05
) {
    (abs(get_xr(plt_dims)[1]) * portion) - (2 * rad_decrease)
}

# given the circle placements, estimate the legend dataframe with the least
# *number* of circles that overlap
estimate_best_legend_df <- function(
    plt_dims, plot_df, unpositioned_legend_df, legend_dims, buffer
) {
    min_num_circles_covered <- Inf
    best_legend_df <- data.frame()

    for (pos in c("top_left", "top_right", "bottom_left", "bottom_right")) {
        
        curr_legend_df <- generate_legend_df(
            pos, unpositioned_legend_df, legend_dims, plt_dims, buffer
        )

        # TODO this should be estimate from the plot dataframe instead
        curr_num_circles_covered <- num_circles_covered_by_legend(
            plot_df,
            minmax_dims = get_legend_backing_minmax_dims(curr_legend_df)
        )

        if (curr_num_circles_covered == 0) {
            return(curr_legend_df)
        }

        if (curr_num_circles_covered < min_num_circles_covered) {
            min_num_circles_covered <- curr_num_circles_covered
            best_legend_df <- curr_legend_df
        }

    }

    best_legend_df
}

num_circles_covered_by_legend <- function(df, minmax_dims) {
    num_circles_covered <- 0
    for (i in seq_len(nrow(df))) {
        does_overlap_x <- is_bound_between(
            df$x[i],
            lowerbound = minmax_dims["xmin"] + df$r[i],
            upperbound = minmax_dims["xmax"] - df$r[i]
        )

        does_overlap_y <- is_bound_between(
            df$y[i],
            lowerbound = minmax_dims["ymin"] + df$r[i],
            upperbound = minmax_dims["ymax"] - df$r[i]
        )

        num_circles_covered <- num_circles_covered +
            (does_overlap_x && does_overlap_y)
    }

    num_circles_covered
}

gen_unpositioned_legend_df <- function(
    legend_sizes, spacing, circ_scale_factor, rad_decrease, color,
    additional_middle_spacing
) {
    radii <- (sqrt(legend_sizes) * circ_scale_factor) - rad_decrease
    num_radii <- length(radii)
    label_x <- spacing + 1 - rad_decrease +
        (sqrt(radii[num_radii]) * circ_scale_factor) +
        additional_middle_spacing

    data.frame(
        "circle_x" = rep(0, num_radii),
        "label_x" = rep(label_x, num_radii),
        "y" = get_unpositioned_y_coords(radii, spacing),
        "rad" = radii,
        "labels" = as.character(legend_sizes),
        "color" = rep(color, num_radii)
    )
}

get_unpositioned_y_coords <- function(radii, spacing) {
    y_coords <- numeric(length(radii))
    curr_y <- -spacing
    r <- 0

    for (i in seq_along(radii)) {
        prev_radius <- r
        r <- radii[i]
        curr_y <- curr_y - prev_radius - r - spacing
        y_coords[i] <- curr_y
    }

    y_coords
}

# get the width and height of the circles and labels, NOT accounting for buffer
get_legend_dims <- function(unpositioned_legend_df) {
    c(
        "x" = max_rad(unpositioned_legend_df) +
            get_label_x(unpositioned_legend_df),

        "y" = min_rad(unpositioned_legend_df) +
            abs(max_y(unpositioned_legend_df)) -
            abs(min_y(unpositioned_legend_df)) +
            max_rad(unpositioned_legend_df)
    )
}

generate_legend_df <- function(
    pos, unpositioned_legend_df, legend_dims, plt_dims, buffer
) {
    if (!is.numeric(pos)) {
        pos <- estimate_top_left_circ_coord(
            unpositioned_legend_df = unpositioned_legend_df,
            legend_dims = legend_dims,
            destination_str = pos,
            plt_dims = plt_dims,
            buffer = buffer
        )
    }

    move_unpositioned_legend_df(
        unpositioned_legend_df,
        to_top_left_destination_coord = pos
    )
}

min_y <- function(legend_df) legend_df[1, "y"]
max_y <- function(legend_df) legend_df[nrow(legend_df), "y"]

min_rad <- function(legend_df) legend_df[1, "rad"]
max_rad <- function(legend_df) legend_df[nrow(legend_df), "rad"]

get_circle_x <- function(legend_df) legend_df[1, "circle_x"]
get_label_x <- function(legend_df) legend_df[1, "label_x"]

# get the starting coordinate for the top left center of the *first circle*
estimate_top_left_circ_coord <- function(
    unpositioned_legend_df, legend_dims, destination_str, plt_dims, buffer
) {

    xr <- get_xr(plt_dims)
    yr <- get_yr(plt_dims)
    
    if (identical(destination_str, "top_left")) {
        return(c(xr[1] + max_rad(unpositioned_legend_df) + buffer, yr[2]))
    }

    if (identical(destination_str, "top_right")) {
        return(c(xr[2] - unpositioned_legend_df[1, "label_x"] - buffer, yr[2]))
    }

    if (identical(destination_str, "bottom_left")) {
        return(c(
            xr[1] + max_rad(unpositioned_legend_df) + buffer,
            yr[1] + legend_dims[2]
        ))
    }
    # bottom right
    return(c(
        xr[2] - unpositioned_legend_df[1, "label_x"] - buffer,
        yr[1] + legend_dims[2]
    ))
}

move_unpositioned_legend_df <- function(
    unpositioned_df, to_top_left_destination_coord
) {
    dx <- to_top_left_destination_coord[1] - unpositioned_df[1, "circle_x"]
    dy <- to_top_left_destination_coord[2] - unpositioned_df[1, "y"]

    unpositioned_df[, "circle_x"] <- unpositioned_df[, "circle_x"] + dx
    unpositioned_df[, "label_x"] <- unpositioned_df[, "label_x"] + dx
    unpositioned_df[, "y"] <- unpositioned_df[, "y"] + dy

    unpositioned_df
}

# get the coordinate in the middle of the top of the legend
get_legend_title_coord <- function(legend_df, legend_dims, spacing) {
    c("x" = sum(get_legend_backing_minmax_dims(legend_df)[1:2]) * 0.5,
      "y" = min_y(legend_df) + min_rad(legend_df) + (1.5 * spacing))
}

add_legend_backing <- function(plt, legend_df) {

    dims <- get_legend_backing_minmax_dims(legend_df)

    plt + ggplot2::annotate(
        geom = "rect",
        xmin = dims["xmin"], xmax = dims["xmax"],
        ymin = dims["ymin"], ymax = dims["ymax"],
        fill = "white",
        colour = "black"
    )

}

get_legend_backing_minmax_dims <- function(legend_df) {
    spacing <- 0.15
    max_radius <- max_rad(legend_df)
    
    xmin <- get_circle_x(legend_df) - max_radius - spacing
    xmax <- get_label_x(legend_df) + max_radius + spacing

    ymin <- max_y(legend_df) - max_radius - spacing
    ymax <- min_y(legend_df) + min_rad(legend_df) + spacing

    c("xmin" = xmin, "xmax" = xmax, "ymin" = ymin, "ymax" = ymax)
}

# could put the ggplot color legend by sticking some points under something
