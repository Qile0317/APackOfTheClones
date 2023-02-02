#The API functions for generating the ggplot, with imperfect repulsion

#full join a list of lists of (x,y,r) vectors into a dataframe with generated labels. memory inefficient
df_full_join <- function(clstr_list) {
  df <- data.frame(label = character(0),
                   x = numeric(0),
                   y = numeric(0),
                   r = numeric(0))

  for(i in 1:length(clstr_list)){
    if (!is.null(clstr_list[[i]])) {

      df <- invisible(full_join(
        df, data.frame(label = paste("cluster", as.character(i-1)),
                       x = clstr_list[[i]]$x,
                       y = clstr_list[[i]]$y,
                       r=clstr_list[[i]]$rad)))
    }

  }
  return(df)
}

#API for plotting size vectors . Integrating the data is work in progress.
plot_API <- function(sizes, # list of size vectors,[[1]] c(a,b,..)
                     centroids, # centroids of size vectors [[1]] c(x,y)
                     ORDER = TRUE,
                     try_place = TRUE,
                     progbar = TRUE, # packing
                     repulse = FALSE,
                     thr = 1, G = 0.05, max_repulsion_iter = 100, #repulsion parameters
                     n = 360, linetype = "blank",
                     plot_title = "Sizes of clones within each cluster",
                     haslegend=TRUE, void=TRUE, origin=FALSE){ #prev 3 lines were for plotting function
  ans <- list()

  #circle layout
  for(i in 1:length(sizes)){
    if (length(sizes[[i]]) == 0) {
      ans[[i]] <- NULL
    }else{
      if(progbar){print(paste("packing cluster", as.character(i)))}
      ans[[i]] <- circle_layout(sizes[[i]],
                                centroid = centroids[[i]],
                                ORDER = ORDER,
                                try_place = try_place,
                                progbar = progbar)
    }
  }

  #repulsion, not sure how it handles nulls
  if (repulse) {
    if(progbar){print("repulsing clusters")}
    ans <- repulse_cluster(ans, thr = thr, G = G, max_iter = max_repulsion_iter)
  }

  print(head(ans)) # debug

  #joining list into df
  ans <- df_full_join(ans)

  #return(ans) # debug. but everything works untill this point

  #plotting
  ans <- plot_clusters(ans, n = n, linetype = linetype, title = plot_title,
                       haslegend = haslegend, void = void, origin = origin)
  return(ans)
}

#the following is unfinished. this will be the main function that the package was built for and the other two functions can go when this is done
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
                       plot_title = "Sizes of clones within each cluster",
                       has_void_theme = TRUE,
                       has_legend = TRUE,
                       show_origin = FALSE,
                       clone_scale_factor = 0.001) {

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

  # pack and return the plot
  return(plot_API(sizes = clone_size_list,
                  centroids = centroid_list,
                  ORDER = ORDER,
                  try_place = try_place,
                  progbar = progbar,
                  repulse = repulse,
                  thr = repulsion_threshold,
                  G = repulsion_strength,
                  max_repulsion_iter = max_repulsion_iter,
                  n = res,
                  plot_title = plot_title,
                  haslegend = has_legend,
                  void = has_void_theme,
                  origin = show_origin))
}
