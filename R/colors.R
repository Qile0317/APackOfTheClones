# wrapper to get the number of identified clusters:
count_umap_clusters <- function(seurat_obj) {
  return(length(levels(seurat_obj@meta.data[["seurat_clusters"]])))
}

# get the ggplot colors - important function for 'apotc'
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# generate a hashmap of cluster to color
gen_cluster_color_hashmap <- function(num_clusters) {
  color_vec <- gg_color_hue(num_clusters)
  output <- hash::hash()
  for (i in 1:num_clusters) {
    cluster_str <- paste("cluster", as.character(i-1))
    output[cluster_str] <- color_vec[i]
  }
  output
}

#' inserts a list of colors into a column in the cluster df by label with the
#' v0.1.2 version of the packing algos
#' @noRd
#' @importFrom dplyr %>% 
insert_colors <- function(cluster_dataframe, num_clusters) {
  color_hashmap <- gen_cluster_color_hashmap(num_clusters)
  color_vec <- cluster_dataframe[[1]] # 1 is label
  for (i in 1:length(color_vec)) {
    color_vec[i] <- color_hashmap[[color_vec[i]]]
  }
  return(cluster_dataframe %>% dplyr::mutate("color" = color_vec))
}

# in the future should probably make a fake ggplot cluster legend on the right
# side by inserting scatterplt points? (like in the seurat UMPA plot) But also, 
# the problem is that it becomes inconsistent with the clone size legend :/