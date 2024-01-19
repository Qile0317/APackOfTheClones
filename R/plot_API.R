create_initial_apotc_plot <- function(apotc_obj, res, linetype) {
    plot_clusters(
        clusters = get_plottable_df_with_color(apotc_obj),
        n = res,
        linetype = linetype
    )
}

get_plottable_df_with_color <- function(apotc_data) {
    extract_and_add_colors(
      apotc_obj = apotc_data,
      plot_df = df_full_join(get_clusterlists(apotc_data))
    )
}

# full join a list of lists of (x,y,r) vectors into a dataframe with
# generated labels.
df_full_join <- function(clstr_list) {
    df <- data.frame(
        'label' = character(0),
        'x' = numeric(0),
        'y' = numeric(0),
        'r' = numeric(0)
    )

    seurat_cluster_index <- 0 # zero indexed :P
    for (i in seq_along(clstr_list)) {
        if (!isnt_empty_nor_na(clstr_list[[i]])) {
            seurat_cluster_index <- seurat_cluster_index + 1
            next
        }
        df <- dplyr::full_join(
            df,
            data.frame(
                "label" = rep( # TODO this could be customized or add new col
                  paste("cluster", seurat_cluster_index),
                  length(clstr_list[[i]][["x"]])
                ),
                "x" = clstr_list[[i]][["x"]],
                "y" = clstr_list[[i]][["y"]],
                "r" = clstr_list[[i]][["rad"]]
            ),
            by = dplyr::join_by("label", "x", "y", "r")
        )
        seurat_cluster_index <- seurat_cluster_index + 1
    }
    df
}

# result plotting function. clusters is a list of clusterlists TRANSFORM into a
# dataframe, which are clusters. A cluster list includes x, y, rad, centroid,
# clRad. the clusters imput is a dataframe.

plot_clusters <- function(
  clusters,
  n = 360,
  linetype = "blank"#,
  #alpha = 1
) {
  ggplot2::ggplot(data = clusters) +
    ggforce::geom_circle(
      mapping = apotc_aes_string(
        x0 = "x",
        y0 = "y",
        r = "r",
        fill = "color"#,
        #alpha = as.character(alpha)
      ),
      n = n,
      linetype = linetype
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme(legend.position = "none")
}

# TODO not quite the same
add_default_theme <- function(plt, reduction) {
	label_hashmap <- hash::hash(
		c("umap", "tsne", "pca"), c("UMAP", "tSNE", "PC")
	)
	label <- label_hashmap[[reduction]]

	plt +
		ggplot2::theme_classic() +
		ggplot2::xlab(paste(label, 1, sep = "_")) +
		ggplot2::ylab(paste(label, 2, sep = "_"))
}

# TODO make better - should probably scale instead
get_retain_scale_dims <- function(
  seurat_obj, reduction, ball_pack_plt, plot_dims
) {

  reduction_dims <- get_plot_dims(
    Seurat::DimPlot(object = seurat_obj, reduction = reduction)
  )
  reduction_xr <- get_xr(reduction_dims)
  reduction_yr <- get_yr(reduction_dims)

  ball_pack_xr <- get_xr(plot_dims)
  ball_pack_yr <- get_yr(plot_dims)

  # set new ranges
  min_xr <- min(ball_pack_xr[1], reduction_xr[1])
  max_xr <- max(ball_pack_xr[2], reduction_xr[2])

  min_yr <- min(ball_pack_yr[1], reduction_yr[1])
  max_yr <- max(ball_pack_yr[2], reduction_yr[2])

  # return dims in same output format as get_plot_dims
  list("xr" = c(min_xr, max_xr), "yr" = c(min_yr, max_yr))
  
}
