# function to mutate the seurat object `apotc` reduction attribute
# and incorporate the clusters, instead of returning



# in the future, make it also possible with PCA and TSNE, not to mention
# multi-sample data, BCR, and SCE objects from scRepertoire.

pack_clonal_clusters <- function(
  integrated_seurat_obj,
  rad_decrease = 1,
  ORDER = TRUE,
  try_place = TRUE,
  progbar = TRUE){
  
  ans <- list()
  clonedat <- integrated_seurat_obj@reductions[["apotc"]]
  sizes <- sqrt(clonedat@clone_sizes) * clonedat@clone_scale_factor

  #circle layout
  for(i in 1:clonedat@num_clusters){
    if (is.null(sizes[[i]])) {
      ans[[i]] <- NA # important!
    }else{
      if(progbar){
        message(paste("\npacking cluster", as.character(i-1)))
      }
      ans[[i]] <- circle_layout(
        sizes[[i]],
        centroid = clonedat@centroids[[i]],
        rad_decrease = rad_decrease,
        ORDER = ORDER,
        try_place = try_place,
        progbar = progbar
      )
    }
  }
  integrated_seurat_obj
}
