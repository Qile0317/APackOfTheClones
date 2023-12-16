# script to make a custom circle size legend overlay
# needs a rehaul to use an ApotcData getter to get the sizes incase the
# implemntation changes in the future (basic OOP principle :P)

estimate_legend_sizes <- function(apotc_obj) {
    sizes <- unlist(as.numeric(apotc_obj@clone_sizes))
    sort(unique(round(c(1, median(sizes), mean(sizes), max(sizes)))))
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
    do_add_legend_border = FALSE # TODO
) {

    # setup relevant variables
    rad_decrease <- get_rad_decrease(apotc_obj)

    if (should_estimate(spacing))
        spacing <- calculate_legend_spacing(spacing, plt, rad_decrease)

    if (should_estimate(sizes))
        sizes <- estimate_legend_sizes(apotc_obj)
    
    if (!is.numeric(pos))
        pos <- correct_legend_coord_str(pos)

    # calculate circle positions

    unpositioned_legend_df <- gen_unpositioned_legend_df(
        legend_sizes = sizes,
        spacing = spacing,
        circ_scale_factor = apotc_obj@clone_scale_factor,
        rad_decrease = rad_decrease,
        color = color
    )

    legend_dims <- get_legend_dims(unpositioned_legend_df)

    if (!is.numeric(pos)) {
        pos <- estimate_top_left_circ_coord(
            unpositioned_legend_df = unpositioned_legend_df,
            legend_dims = legend_dims,
            destination_str = pos,
            plt = plt,
            buffer = buffer
        )
    }

    legend_df <- move_unpositioned_legend_df(
        unpositioned_legend_df,
        to_top_left_destination_coord = pos
    )

    label_coord <- get_legend_title_coord(legend_df, spacing)

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

correct_legend_coord_str <- function(pos) {
    pos <- tolower(strip_spaces(pos))
    pos <- switch(
        pos,
        "topleft" = "top_left",
        "topright" = "top_right",
        "bottomleft" = "bottom_left",
        "bottomright" = "bottom_right",
        pos
    )

    user_attempt_correction(
        s = pos,
        strset = c("top_left", "top_right", "bottom_left", "bottom_right"),
        stop_msg_start = "invalid legend coordinate string"
    )
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
    y_coords <- vector("numeric", length(radii))
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

get_legend_dims <- function(unpositioned_legend_df) {
    num_circles <- nrow(unpositioned_legend_df)
    
    c(
        "x" = max_rad(unpositioned_legend_df) +
            unpositioned_legend_df[num_circles, "label_x"],

        "y" = unpositioned_legend_df[1, "rad"] +
            unpositioned_legend_df[num_circles, "y"] +
            max_rad(unpositioned_legend_df)
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

    destination_str <- correct_legend_coord_str(destination_str)

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
            yr[1] - legend_dims[2]
        ))
    }
    # bottom right
    return(c(
        xr[2] - unpositioned_legend_df[1, "label_x"] - buffer,
        yr[1] - legend_dims[2]
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
get_legend_title_coord <- function(legend_df, spacing) {
    c(
        "x" = abs(
            legend_df[1, "label_x"] -
                legend_df[1, "circle_x"] - max_rad(legend_df)
        ),
        "y" = legend_df[1, "y"] - spacing
    )
}

add_legend_backing <- function(plt, legend_df) {
    spacing <- 0.15
    linewidth <- bound_num(
        abs(get_xr(plt)[2] - get_xr(plt)[1]) * 0.002, 0.001, 0.1
    )
    max_radius <- max_rad(legend_df)

    xmin <- get_circle_x(legend_df) - max_radius - spacing
    xmax <- get_label_x(legend_df) + max_radius + spacing

    ymin <- max_y(legend_df) + max_radius + spacing
    ymax <- min_y(legend_df) - min_rad(legend_df) - spacing

    # add the back border rectangle
    plt <- plt + ggplot2::geom_rect(aes(
            xmin = xmin - linewidth, xmax = xmax + linewidth,
            ymin = ymin + linewidth, ymax = ymax - linewidth,
            fill = "black"
        ))
    
    # add the white inside
    plt + ggplot2::geom_rect(aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
            fill = "white"
        ))
}

# could put the ggplot color legend by sticking some points under something
