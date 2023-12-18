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

get_plottable_df_with_color <- function(apotc_data) {
    extract_and_add_colors(apotc_data, df_full_join(apotc_data@clusters))
}

# result plotting function. clusters is a list of clusterlists TRANSFORM into a
# dataframe, which are clusters. A cluster list includes x, y, rad, centroid,
# clRad. the clusters imput is a dataframe.

plot_clusters <- function(
  clusters, n = 360, linetype ="blank",
  title = "Sizes of clones within each cluster",
  haslegend = FALSE, void = TRUE,
  origin = FALSE
){

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
  } else {
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

add_default_theme <- function(plt, reduction) {
	label_hashmap <- hash::hash(
		c("umap", "tsne", "pca"), c("UMAP", "tSNE", "PC")
	)
	label <- label_hashmap[[reduction]]

	plt +
		ggplot2::theme_classic() +
		ggplot2::xlab(paste(label, 1, sep = "_")) +
		ggplot2::ylab(paste(label, 2, sep = "_")) +
		ggplot2::ggtitle("Sizes of clones within each cluster")
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
