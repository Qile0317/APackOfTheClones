# script to make a custom circle size legend overlay
# TODO estimate optimal legend placement position

estimate_legend_sizes <- function(apotc_obj) {
    sizes <- unlist(apotc_obj@clone_sizes)
    sort(unique(round(c(
        1, median(sizes), mean(sizes), mean(unique(sizes)), max(sizes)
    ))))
}

insert_legend <- function(
    plt,
    apotc_obj,
    sizes,
    pos,
    buffer,
    color = "#808080",
    n = 360,
    spacing = "auto",
    legend_label = "Clone sizes",
    legend_textsize = 5,
    do_add_legend_border = FALSE
) {

    # setup relevant variables
    rad_decrease <- get_rad_decrease(apotc_obj)

    if (should_estimate(spacing))
        spacing <- calculate_legend_spacing(spacing, plt, rad_decrease)

    if (should_estimate(sizes))
        sizes <- estimate_legend_sizes(apotc_obj)
    
    pos <- correct_legend_coord_if_str(pos)

    # calculate relevant legend plotting data

    unpositioned_legend_df <- gen_unpositioned_legend_df(
        legend_sizes = sizes,
        spacing = spacing,
        circ_scale_factor = apotc_obj@clone_scale_factor,
        rad_decrease = rad_decrease,
        color = color
    )

    legend_dims <- get_legend_dims(unpositioned_legend_df)

    if (should_estimate(pos)) {
        legend_df <- estimate_best_legend_df(
            plt, apotc_obj, unpositioned_legend_df, legend_dims, buffer
        )
    } else {
        legend_df <- generate_legend_df(
            pos, unpositioned_legend_df, legend_dims, plt, buffer
        )
    }

    label_coord <- get_legend_title_coord(legend_df, legend_dims, spacing)

    # plotting

    # add the legend label on top
    plt <- plt + ggplot2::annotate(
        "text", x = label_coord[1], y = label_coord[2],
        label = legend_label, size = legend_textsize
    )

    # add the background
    if (do_add_legend_border) {
        plt <- add_legend_backing(plt, legend_df)
    }
    
    # add the side number labels
    plt <- plt + ggplot2::annotate(
        "text", x = legend_df[, "label_x"], y = legend_df[, "y"],
        label = legend_df[, "labels"], size = legend_textsize
    )

    # add the circles and return
    plt + ggforce::geom_circle(
        data = legend_df,
        mapping = apotc_aes_string(
            x0 = "circle_x",
            y0 = "y",
            r = "rad",
            fill = "color"
        ),
        linetype = "blank",
        n = n
    )
}

calculate_legend_spacing <- function(
    spacing, plt, rad_decrease, portion = 0.05
) {
    (abs(get_xr(plt)[1]) * portion) - (2 * rad_decrease)
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

# given the circle placements, estimate the legend dataframe with the least
# *number* of circles that overlap
estimate_best_legend_df <- function(
    plt, apotc_obj, unpositioned_legend_df, legend_dims, buffer
) {
    min_num_circles_covered <- Inf
    best_legend_df <- data.frame()

    for (pos in c("top_left", "top_right", "bottom_left", "bottom_right")) {

        curr_legend_df <- generate_legend_df(
            pos, unpositioned_legend_df, legend_dims, plt, buffer
        )

        minmax_dims <- get_legend_backing_minmax_dims(plt, curr_legend_df)

        curr_num_circles_covered <- num_circles_covered_by_legend(
            apotc_obj, minmax_dims
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

num_circles_covered_by_legend <- function(apotc_obj, minmax_dims) {
    num_circles_covered <- 0
    for (cluster in apotc_obj@clusters) {
        for (i in seq_along(cluster$x)) {
            does_overlap_x <- is_bound_between(
                cluster$x[i],
                lowerbound = minmax_dims["xmin"] + cluster$rad[i],
                upperbound = minmax_dims["xmax"] - cluster$rad[i]
            )

            does_overlap_y <- is_bound_between(
                cluster$y[i],
                lowerbound = minmax_dims["ymin"] + cluster$rad[i],
                upperbound = minmax_dims["ymax"] - cluster$rad[i]
            )

            if (does_overlap_x && does_overlap_y) {
                num_circles_covered <- num_circles_covered + 1
            }
        }
        
    }
    num_circles_covered
}

gen_unpositioned_legend_df <- function(
    legend_sizes, spacing, circ_scale_factor, rad_decrease, color
) {
    radii <- (sqrt(legend_sizes) * circ_scale_factor) - rad_decrease
    num_radii <- length(radii)
    label_x <- spacing + 1 - rad_decrease +
        (sqrt(radii[num_radii]) * circ_scale_factor)

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
    pos, unpositioned_legend_df, legend_dims, plt, buffer
) {
    if (!is.numeric(pos)) {
        pos <- estimate_top_left_circ_coord(
            unpositioned_legend_df = unpositioned_legend_df,
            legend_dims = legend_dims,
            destination_str = pos,
            plt = plt,
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
    unpositioned_legend_df, legend_dims, destination_str, plt, buffer
) {

    xr <- get_xr(plt)
    yr <- get_yr(plt)
    
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
    c("x" = get_circle_x(legend_df) - max_rad(legend_df) + (legend_dims[1] / 2),
      "y" = min_y(legend_df) + min_rad(legend_df) + spacing)
}

add_legend_backing <- function(plt, legend_df) {
    linewidth <- get_linewidth(plt)
    dims <- get_legend_backing_minmax_dims(plt, legend_df)

    # add the back border rectangle
    plt <- plt + ggplot2::geom_rect(aes(
            xmin = dims["xmin"] - linewidth, xmax = dims["xmax"] + linewidth,
            ymin = dims["ymin"] - linewidth, ymax = dims["ymax"] + linewidth,
            fill = "black"
        ))
    
    # add the white inside
    plt + ggplot2::geom_rect(aes(
            xmin = dims["xmin"], xmax = dims["xmax"],
            ymin = dims["ymin"], ymax = dims["ymax"],
            fill = "white",
            linetype = "blank"
        )) +
        ggplot2::theme(legend.position = "none")
}

get_linewidth <- function(plt) {
    xr <- get_xr(plt)
    bound_num(
        abs(xr[2] - xr[1]) * 0.002,
        lowerbound = 0.001,
        upperbound = 0.1
    )
}

get_legend_backing_minmax_dims <- function(plt, legend_df) {
    spacing <- 0.15
    max_radius <- max_rad(legend_df)
    
    xmin <- get_circle_x(legend_df) - max_radius - spacing
    xmax <- get_label_x(legend_df) + max_radius + spacing

    ymin <- max_y(legend_df) - max_radius - spacing
    ymax <- min_y(legend_df) + min_rad(legend_df) + spacing

    c("xmin" = xmin, "xmax" = xmax, "ymin" = ymin, "ymax" = ymax)
}

# could put the ggplot color legend by sticking some points under something
