create_initial_apotc_plot <- function(apotc_obj, res, linetype, alpha) {
    plot_clusters(
        clusters = get_plottable_df_with_color(apotc_obj),
        num_total_clusters = get_num_clusters(apotc_obj),
        n = res,
        linetype = linetype,
        alpha = alpha
    )
}

get_plottable_df_with_color <- function(apotc_data) {
    extract_and_add_colors(
        apotc_obj = apotc_data,
        plot_df = df_full_join(get_clusterlists(apotc_data))
    )
}

# full join a list of clusterlists vectors into a dataframe with
# generated labels.
df_full_join <- function(clstr_list) {

    df <- data.frame(
        "label" = character(0),
        "x" = numeric(0),
        "y" = numeric(0),
        "r" = numeric(0)
    ) %>%
        update_clusterlist_df(get_first_nonempty(clstr_list))
    
    if (contains_clonotypes(get_first_nonempty(clstr_list))) {
        cols_to_join_by <- dplyr::join_by("label", "x", "y", "r", "clonotype")
    } else {
        cols_to_join_by <- dplyr::join_by("label", "x", "y", "r")
    }

    seurat_cluster_index <- 0

    for (i in seq_along(clstr_list)) {

        seurat_cluster_index <- seurat_cluster_index + 1
        if (is_empty(clstr_list[[i]])) next

        df <- df %>% dplyr::full_join(
            convert_to_dataframe(clstr_list[[i]], seurat_cluster_index - 1), # zero indexed
            by = cols_to_join_by
        )
    }

    df
}

plot_clusters <- function(
  clusters,
  num_total_clusters,
  n = 360,
  linetype = "blank",
  alpha = 1
) {
  clusters %>%
    ggplot2::ggplot() +
      ggforce::geom_circle(
        apotc_aes_string(
          x0 = "x",
          y0 = "y",
          r = "r",
          fill = "color"
        ),
        n = n,
        linetype = linetype,
        alpha = alpha
      ) +
      ggplot2::coord_fixed() +
      ggplot2::scale_fill_identity()
}

add_default_theme <- function(plt, reduction) {
	label_hashmap <- hash::hash(
		  c("umap", "tsne", "pca"), c("UMAP", "tSNE", "PC")
	)
	label <- label_hashmap[[reduction]]
  if (is.null(label)) label <- reduction

	plt +
      ggplot2::theme_classic() +
      ggplot2::xlab(paste(label, 1, sep = "_")) +
      ggplot2::ylab(paste(label, 2, sep = "_"))
}

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

check_is_apotc_ggplot <- function(x) {
  if (!isApotcGGPlot(x)) {
    stop(call. = FALSE, "not an output of `APOTCPlot` or `vizAPOTC`")
  }
}

remove_ggplot_layers <- function(ggplot_obj, layer_indicies) {
  ggplot_obj$layers[layer_indicies] <- NULL
  ggplot_obj
}

get_ggplot_data <- function(x) x$data

set_ggplot_data <- function(ggplot_obj, new_data) {
  ggplot_obj$data <- new_data
  ggplot_obj
}
