#' @title 
#' Calculate seurat cluster centroids based on a Dimensional reduction
#' 
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Utility function to calculate the physical xy coordinates of each seurat
#' cluster based on a dimensional reduction already present in the object.
#' The results are returned in a list with the length of the number of distinct
#' seurat clusters based on the seurat_obj `meta.data`.
#' 
#' @param seurat_obj input seurat object with the dimensional reduction of
#' choice already present, and seurat clusters computed.
#' @param reduction character. The reduction that the centroid calculation
#' should be based on. Currently, can only be "umap", "tsne", or "pca".
#' 
#' @return A list of the length of the number of distinct clusters in the
#' seurat object metadata, where each element of the list is a numeric vector
#' of length 2, with the numbers corresponding to the x and y coordinate
#' respectively of the seurat cluster with the corresponding index.
#' 
#' @export
#' 
#' @examples
#' data("combined_pbmc")
#' getReductionCentroids(combined_pbmc, reduction = "umap")
#' #>
#' 
getReductionCentroids <- function(seurat_obj, reduction = "umap") {
  get_cluster_centroids(
    seurat_obj = seurat_obj,
    reduction = user_get_reduc_obj(seurat_obj, reduction),
    passed_in_reduc_obj = TRUE
  )
}

# Returns the number of valid barcodes that are not NA's
count_tcr_barcodes <- function(seurat_obj) {
  sum(!is.na(seurat_obj@meta.data[["barcode"]]))
}

count_clones <- function(seurat_obj, clonecall) {
  sum(!is.na(seurat_obj@meta.data[[clonecall]]))
}

# get the percent of NA's in the metadata barcode column for the message
percent_na <- function(seurat_obj) {
  num_barcodes <- length(seurat_obj@meta.data[["barcode"]])
  100 * (num_barcodes - count_tcr_barcodes(seurat_obj)) / num_barcodes
}

get_rna_assay_barcodes <- function(seurat_obj) {
    seurat_obj@assays[["RNA"]]@data@Dimnames[[2]]
}

# seurat cluster related functions

count_num_clusters <- function(seurat_obj) {
  data.table::uniqueN((seurat_obj@meta.data[["seurat_clusters"]]))
}

get_num_total_clusters <- function(seurat_obj) {
  length(levels(seurat_obj@meta.data[["seurat_clusters"]]))
}

# reduction related functions

get_2d_embedding <- function(seurat_obj, reduction) {
  seurat_obj@reductions[[reduction]]@cell.embeddings[, 1:2]
}

attempt_correction <- function(s) {
    user_attempt_correction(
      s = strip_spaces(tolower(s)),
      strset = c("umap", "tsne", "pca"),
      stop_msg_start = "Invalid reduction"
    )
}

is_reduction_name <- function(reduction) {
  any(reduction == c("umap", "tsne", "pca"))
}

user_get_reduc_obj <- function(seurat_obj, reduction) {
  if (!is_seurat_object(seurat_obj)) {
    stop("not a seurat object!")
  }

  reduction <- attempt_correction(reduction)

  if (!is_reduction_name(reduction)) {
    stop(paste("A", reduction, "has not been ran on the seurat object"))
  }

  seurat_obj@reductions[[reduction]]
}
