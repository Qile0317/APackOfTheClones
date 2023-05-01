# script to make a custom circle size legend overlay

library(ggplot2)
library(ggforce)

# place the x0 and y0 of the "Clone sizes" text label. Currently bugged.
# In the future it can be more automated and an optional rectangle can be added
#' @noRd
#' @import ggforce
insert_legend <- function(plt, circ_scale_factor, sizes = c(1,5,10), pos = "top_left", buffer = 0.5, color = "#808080",n=360) {
  
  xr <- ggplot2::ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range
  yr <- ggplot2::ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range
  
  # this is the most rudimentary legend placement; purely on the corner.
  pos_list <- list(
    "top_left" = c(xr[1] + buffer, yr[2] - buffer),
    "top_right" = c(xr[2] - buffer, yr[2] - buffer),
    "bottom_left" = c(xr[1] + buffer, yr[1] + buffer),
    "bottom_right" = c(xr[2] - buffer, yr[1] + buffer)
  )
  coord <- pos_list[[pos]]
  plt <- plt + ggplot2::annotate("text", x = coord[1], y = coord[2], label = "Clone sizes")
  coord[2] <- coord[2] - 0.5
  r <- 0
  for (size in sort(unique(sizes))) {
    prev_r <- r
    r <- sqrt(size)*circ_scale_factor
    coord[2] <- coord[2] - prev_r - r - 0.1
    print(coord)
    plt <- plt + geom_circle(ggplot2::aes(x0=coord[1],y0=coord[2],r=r,fill=color),linetype ="blank",n=n)
    plt <- plt + ggplot2::annotate("text", x = coord[1] + r + 0.75, y = coord[2], label = as.character(size))
  }
  return(plt)
}
