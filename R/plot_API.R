# full join a list of lists of (x,y,r) vectors into a dataframe with
# generated labels.
df_full_join <- function(clstr_list) {
    df <- data.frame(
        'label' = character(0),
        'x' = numeric(0),
        'y' = numeric(0),
        'r' = numeric(0)
    )

    seurat_cluster_index <- 0
    for(i in seq_along(clstr_list)){
        if (isnt_empty_nor_na(clstr_list[[i]])) {
            df <- dplyr::full_join(
                df,
                data.frame(
                  "label" = rep(
                    paste("cluster", seurat_cluster_index),
                    length(clstr_list[[i]][["x"]])
                  ),
                  "x" = clstr_list[[i]][["x"]],
                  "y" = clstr_list[[i]][["y"]],
                  "r" = clstr_list[[i]][["rad"]]
                ),
                by = dplyr::join_by("label", "x", "y", "r")
            )
        }
        seurat_cluster_index <- seurat_cluster_index + 1
    }
    df
}


# result plotting function. clusters is a list of clusterlists TRANSFORM into a
# dataframe, which are clusters. A cluster list includes x, y, rad, centroid,
# clRad. the clusters imput is a dataframe.

plot_clusters <- function(
  clusters, n = 360, linetype ="blank",
  title = "Sizes of clones within each cluster",
  haslegend = FALSE, void = TRUE,
  origin = FALSE){

  if (!origin) {
    p1 <- ggplot2::ggplot(data = clusters) +
      ggforce::geom_circle(
        mapping = apotc_aes_string(
          x0 = "x", y0 = "y", r = "r", fill = "color"
        ),
        n = n,
        linetype = linetype
      ) +
      ggplot2::scale_fill_identity() +
      ggplot2::labs(title = title) +
      ggplot2::coord_fixed()
  }else {
    p1 <- ggplot2::ggplot(
      clusters,
      mapping = apotc_aes_string(x = "x", y = "y")
      ) +
      ggplot2::geom_point() +
      ggplot2::labs(title = title) +
      ggplot2::coord_fixed()
  }

  if (void) {
    p1 <- p1 + ggplot2::theme_void()
  }

  return(p1 + ggplot2::theme(legend.position = "none"))
}

# not rlly needed anymore
plot_API <- function(
  sizes, # list of size vectors,[[1]] c(a,b,..)
  centroids, # centroids of size vectors [[1]] c(x,y)
  num_clusters,
  rad_decrease = 0,
  ORDER = TRUE,
  scramble = FALSE,
  try_place = FALSE,
  progbar = TRUE,
  repulse = FALSE,
  thr = 1, G = 1,
  max_repulsion_iter = 10,
  n = 360, linetype = "blank",
  plot_title = "Sizes of clones within each cluster",
  haslegend = TRUE,
  void = TRUE,
  origin = FALSE,
  debug_mode = FALSE
) {

  ans <- pack_into_clusterlists(
    sizes, centroids, num_clusters, rad_decrease, ORDER, scramble, try_place,
    progbar
  )

  if (repulse) {
    if(progbar){
      message(
        paste("\nrepulsing all clusters | max iterations =", max_repulsion_iter)
      )
    }
    ans <- repulse_cluster(
      ans, thr = thr, G = G, max_iter = max_repulsion_iter, verbose = progbar
    )
  }

  #joining list into df, add color, then plot.
  ans <- df_full_join(ans)
  ans <- insert_colors(ans, num_clusters)

  plt <- plot_clusters(
    ans, n = n, linetype = linetype, title = plot_title, haslegend = haslegend,
    void = void, origin = origin
  )
  plt
}


# change the axis scales to fit the original plot approximately. Looks pretty bad atm.
# A more advanced version could multiply axses by a small amount to retain ratios exactly
retain_scale <- function(seurat_obj, reduction, ball_pack_plt, buffer = 0) {

  test_reduction_plot <- Seurat::DimPlot(
    object = seurat_obj,
    reduction = reduction
  )

  # get current ranges
  reduction_xr <- get_xr(test_reduction_plot)
  reduction_yr <- get_yr(test_reduction_plot)

  rm("test_reduction_plot")

  ball_pack_xr <- get_xr(ball_pack_plt)
  ball_pack_yr <- get_yr(ball_pack_plt)

  # set new ranges
  min_xr <- min(ball_pack_xr[1], reduction_xr[1]) - buffer
  max_xr <- max(ball_pack_xr[2], reduction_xr[2]) + buffer

  min_yr <- min(ball_pack_yr[1], reduction_yr[1]) - buffer
  max_yr <- max(ball_pack_yr[2], reduction_yr[2]) + buffer

  return(
    ball_pack_plt + ggplot2::coord_cartesian(
      xlim = c(min_xr, max_xr),
      ylim = c(min_yr, max_yr)
    )
  )
}

