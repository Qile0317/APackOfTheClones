# script to make a custom circle size legend overlay - need complete rework.
# Should take up a percentage of corner not a constant amount

get_legend_coordinate <- function(plt, pos, buffer) {
    coord <- numeric(2)

    if (is.character(pos)) {
      xr <- get_xr(plt)
      yr <- get_yr(plt)

      pos_list <- list(
          "top_left" = c(xr[1] + buffer, yr[2]),
          "top_right" = c(xr[2] - buffer, yr[2]),
          "bottom_left" = c(xr[1] + buffer, yr[1] + buffer),
          "bottom_right" = c(xr[2] - buffer, yr[1] + buffer)
      )
        coord <- pos_list[[tolower(pos[1])]]
    }else if(is.numeric(pos)) {
        coord <- pos
    }
    coord
}

calculate_legend_spacing <- function(spacing, plt, portion = 0.05) {
    if (should_estimate(spacing)) {
        return(get_xr(plt) * portion)
    }
    spacing
}

insert_legend <- function(
    plt,
    circ_scale_factor,
    rad_decrease,
    sizes = c(1,5,50),
    pos = "top_left",
    buffer = 1.5,
    color = "#808080",
    n = 360,
    spacing = "auto", # spacing should be a percentage of plot height
    legend_label = "Clone sizes",
    legend_textsize = 5
) {

    sizes <- sort(unique(sizes))
    coord <- get_legend_coordinate(plt, pos, buffer)
    spacing <- calculate_legend_spacing(spacing, plt)

    plt <- plt + ggplot2::annotate(
        "text", x=coord[1], y=coord[2], label=legend_label, size=legend_textsize
    )

    # initialize the circles in the legend
    m <- length(sizes)
    legend_list <- list(
        "x" = rep(coord[1], m),
        "y" = c(coord[2] - 0.5, numeric(m-1)),
        "r" = (sqrt(sizes) * circ_scale_factor) - rad_decrease,
        "color" = rep(color, m)
    )

    number_label_x_coord <- coord[1] + spacing + 1 - rad_decrease +
        (sqrt(sizes[m]) * circ_scale_factor)

    coord[2] <- coord[2] - spacing
    r <- 0
    for (i in 1:m) {
        prev_r <- r
        r <- legend_list[[3]][i]
        coord[2] <- coord[2] - prev_r - r - spacing
        legend_list[[2]][i] <- coord[2]

        plt <- plt + ggplot2::annotate(
            "text",
            x = number_label_x_coord,
            y = coord[2],
            label = as.character(sizes[i]),
            size = legend_textsize
        )
    }

    plt + ggforce::geom_circle(
        data = data.frame(legend_list),
        mapping = apotc_aes_string(
            x0 = "x",
            y0 = "y",
            r = "r",
            fill = "color"
        ),
        linetype ="blank",
        n = n
    )
}

# could put the ggplot color legend by sticking some points under something
#insert_color_legend <- function(plt, seurat_obj) {
    #yr <- get_yr(plt)
    #seurat_obj@reductions[["apotc"]]
#}
