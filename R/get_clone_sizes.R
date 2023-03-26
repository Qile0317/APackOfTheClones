# memory and speed inefficient counting of clonotypes within each cluster
get_clone_sizes <- function(integrated_seurat_obj, scale_factor = 0.001) {
  
  df <- data.frame(
    "clusters" = integrated_seurat_obj@meta.data[["seurat_clusters"]],
    "clonotype_ids" = integrated_seurat_obj@meta.data[["raw_clonotype_id"]])
  
  freq <- aggregate(clonotype_ids ~ clusters, data = df, function(x) table(x))
  freq <- freq[["clonotype_ids"]]
  
  for (i in 1:length(freq)) {
    freq[[i]] <- sapply(as.numeric(freq[[i]]), function(x) (sqrt(x) * scale_factor))
  }
  
  return(freq)
}

# ignore the code below it was a prior attemp at implemention an Rcpp version

# wrapper to get the number of identified clusters:
#count_umap_clusters <- function(seurat_obj) {
#  return(length(levels(seurat_obj@meta.data[["seurat_clusters"]])))
#}

#library(Rcpp)
#sourceCpp("src/get_clone_sizes.cpp")

# counting clonotypes to make the sizes of clones within cluster data.
#get_clone_sizes <- function(integrated_seurat_obj, scale_factor = 0.001) {
  #num_clusters <- as.numeric(count_umap_clusters(integrated_seurat_obj))

  ##work in progress
  ##label_radii_index_dict <- dict()
  ##label_list <- list()
  ##color_list <- list()

  ##variables for iteration
  #barcodes <- integrated_seurat_obj@meta.data[["barcode"]]
  #clusters <- as.numeric(integrated_seurat_obj@meta.data[["seurat_clusters"]])
  #clonotype_ids <- integrated_seurat_obj@meta.data[["raw_clonotype_id"]]

  #return(get_clone_sizes_Cpp(
  #  barcodes, clusters, clonotype_ids,
  #  num_clusters, scale_factor
  #  )
  #)
#}

# it can return a list of tabled values for memory efficiency but then some stuff later have to be re-done
