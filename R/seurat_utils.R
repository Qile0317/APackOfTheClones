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

# get the number of seurat clusters
get_num_clusters <- function(seurat_obj) {
  length(levels(seurat_obj@meta.data[["seurat_clusters"]]))
}
