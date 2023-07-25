# script to make a custom circle size legend overlay - need complete rework.
# Should take up a percentage of corner not a constant amount

# inert a circle size legend into a clonal expansion plot
insert_legend <- function(
  plt,
  circ_scale_factor,
  rad_scale_factor,
  sizes = c(1,5,50),
  pos = "top_left",
  buffer = 1.5,
  color = "#808080",
  n = 360,
  spacing = 0.25
) {

  sizes <- sort(unique(sizes))
  xr <- ggplot2::ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range
  yr <- ggplot2::ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range

  if (is.character(pos)) {
    pos_list <- list(
      "top_left" = c(xr[1] + buffer, yr[2]),
      "top_right" = c(xr[2] - buffer, yr[2]),
      "bottom_left" = c(xr[1] + buffer, yr[1] + buffer),
      "bottom_right" = c(xr[2] - buffer, yr[1] + buffer)
    )
    coord <- pos_list[[tolower(pos)]]
  }else if(is.numeric(pos)) {
    coord <- pos
  }

  plt <- plt + ggplot2::annotate(
    "text", x = coord[1], y = coord[2], label = "Clone sizes"
  )

  m <- length(sizes)
  legend_list <- list(
    "x" = rep(coord[1], m),
    "y" = c(coord[2] - 0.5, numeric(m-1)),
    "r" = sqrt(sizes) * circ_scale_factor * rad_scale_factor,
    "color" = rep(color, m)
  )
  coord[2] <- coord[2] - 0.5
  r <- 0
  for (i in 1:m) {
    prev_r <- r
    r <- legend_list[[3]][i]
    coord[2] <- coord[2] - prev_r - r - spacing
    legend_list[[2]][i] <- coord[2]
    plt <- plt + ggplot2::annotate(
      "text",
      x = coord[1] + r + 1,
      y = coord[2],
      label = as.character(sizes[i])
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
    n = n)
}
