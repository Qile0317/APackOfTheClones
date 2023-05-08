# memory and speed inefficient counting of clonotypes within each cluster
# there needs to be testcases for when clusters have no clones
get_clone_sizes <- function(integrated_seurat_obj, scale_factor = 0.001) {
  df <- data.frame(
    "clusters" = integrated_seurat_obj@meta.data[["seurat_clusters"]],
    "clonotype_ids" = integrated_seurat_obj@meta.data[["raw_clonotype_id"]]
  )
  
  freq_df <- aggregate(clonotype_ids ~ clusters, data = df, function(x) table(x))
  
  num_clusters <- length(levels(integrated_seurat_obj@meta.data[["seurat_clusters"]]))
  freq <- vector("list", num_clusters) # each element initialzes to NULL
  for (i in 1:nrow(freq_df)) {
    freq[[freq_df[[1]][i]]] <- sqrt(as.numeric(freq_df[[2]][[i]])) * scale_factor
  }
  freq
}

#' count the number of clonotype sizes per cell cluster in a seurat object integrated with a TCR library
#' 
#' @param integrated_seurat_obj Seurat object that has been integrated with a T-cell receptor library with \code{\link{integrate_tcr}}. More specifically, in the metadata, there must at least be the elements `seurat_clusters` and `raw_clonotype_id`
#' 
#' @return Returns a list of `table` objects, where each element is tabled clonotype frequencies for the seurat cluster corresponding to the same index - 1. For example, the 5th element is a tabled frequency of counts that corresponds to the 4th seurat cluster (as seurat clusters are 0-indexed). If an element is `NULL`, it indicates that there were no corresponding T-cell receptor barcode for the cells in the cluster.
#' 
#' @seealso \code{\link{integrate_tcr}}
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' library(Seurat)
#' library(APackOfTheClones)
#' data("mini_clonotype_data","mini_seurat_obj")
#'
#' # produce an integrated seurat_object
#' integrated_seurat_object <- integrate_tcr(mini_seurat_obj, mini_clonotype_data)
#' clonotype_counts <- count_clone_sizes(integrated_seurat_object)
#' }
count_clone_sizes <- function(integrated_seurat_obj) {
  if (is.null(integrated_seurat_obj@meta.data[["seurat_clusters"]])) {
    stop("A UMAP must first be run on the seurat object")
  }
  if (is.null(integrated_seurat_obj@meta.data[["raw_clonotype_id"]])) {
    stop("Seurat object is not integrated with a T-cell receptor library or has no metadata `raw_clonotype_id`")
  }
  
  df <- data.frame(
    "clusters" = integrated_seurat_obj@meta.data[["seurat_clusters"]],
    "clonotype_ids" = integrated_seurat_obj@meta.data[["raw_clonotype_id"]]
  )
  
  freq_df <- aggregate(clonotype_ids ~ clusters, data = df, function(x) table(x))
  
  num_clusters <- length(levels(integrated_seurat_obj@meta.data[["seurat_clusters"]]))
  freq <- vector("list", num_clusters) # each element initialzes to NULL
  for (i in 1:nrow(freq_df)) {
    freq[[freq_df[[1]][i]]] <- table(as.numeric(freq_df[[2]][[i]]))
  }
  freq
}