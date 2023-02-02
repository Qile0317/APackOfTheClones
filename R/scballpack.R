#' @export
scballpack <- function(seurat_obj, tcr_df,
                       cluster_labels = c(), # will introduce later
                       res = 360,
                       ORDER = TRUE,
                       try_place = TRUE,
                       progbar = TRUE,
                       repulse = TRUE,
                       repulsion_threshold = 1,
                       repulsion_strength = 0.05,
                       max_repulsion_iter = 100,
                       use_default_theme = TRUE,
                       show_origin = FALSE,
                       clone_scale_factor = 0.01,
                       retain_axis_scales = TRUE) {

  # errors/warnings:
  if (is.null(seurat_obj@reductions[["umap"]])) {stop("No UMAP reduction found on the seurat object")}
  if (max_repulsion_iter > 1000) {warning("Repulsion iteration count is high, consider reducing max_repulsion_iter if runtime is too long")}

  # integrate TCR and print how many were integrated
  integrated_seurat_obj <- integrate_tcr(seurat_obj, tcr_df)
  percent_integrated <- 100 - percent_na(integrated_seurat_obj)
  message("")
  message(paste("Percent of cells integrated:", as.character(percent_integrated), "%"))

  # get the clone sizes and cluster centroids
  clone_size_list <- get_clone_sizes(integrated_seurat_obj, scale_factor = clone_scale_factor) # need to scale
  centroid_list <- get_cluster_centroids(integrated_seurat_obj)

  # pack the plot
  result_plot <- plot_API(sizes = clone_size_list,
                          centroids = centroid_list,
                          ORDER = ORDER,
                          try_place = try_place,
                          progbar = progbar,
                          repulse = repulse,
                          thr = repulsion_threshold,
                          G = repulsion_strength,
                          max_repulsion_iter = max_repulsion_iter,
                          n = res,
                          origin = show_origin)

  # retain axis scales on the resulting plot. The function is broken and needs to be fixed
  if (retain_axis_scales) {
    result_plot <- retain_scale(seurat_obj, result_plot)
  }

  #set theme
  if (use_default_theme) {
    result_plot <- result_plot +
      theme_classic() +
      xlab("UMAP_1") +
      ylab("UMAP_2") +
      ggtitle("Sizes of clones within each cluster")
  }

  message("Packing completed successfully. It's highly recommended to use the svglite package to save the resulting plot as an svg")
  return(result_plot)
}
