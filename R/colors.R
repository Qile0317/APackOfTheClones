# wrapper to get the number of identified clusters:
count_umap_clusters <- function(seurat_obj) {
    length(levels(seurat_obj@meta.data[["seurat_clusters"]]))
}

# get the ggplot colors - important function for 'apotc'
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

# generate a hashmap of cluster to color
gen_cluster_color_hashmap <- function(num_clusters) {
  color_vec <- gg_color_hue(num_clusters)
  output <- hash::hash()
  for (i in 1:num_clusters) {
    cluster_str <- paste("cluster", as.character(i - 1))
    output[cluster_str] <- color_vec[i]
  }
  output
}

#' inserts a list of colors into a column in the cluster df by label with the
#' v0.1.2 version of the packing algos
#' @noRd
insert_colors <- function(cluster_dataframe, num_clusters) {
  color_hashmap <- gen_cluster_color_hashmap(num_clusters)
  color_vec <- cluster_dataframe[[1]] # 1 is label
  for (i in seq_along(color_vec)) {
    color_vec[i] <- color_hashmap[[color_vec[i]]]
  }
  return(cluster_dataframe %>% dplyr::mutate("color" = color_vec))
}

# new version that simply takes the readily existing colors in an seurat
# object, and adds the colors to the dataframe

# pair colors to hashmap
pair_colors_to_hash <- function(apotc_obj) {
  color_vec <- apotc_obj@cluster_colors
  output <- hash::hash()
  for (i in 1:apotc_obj@num_clusters) {
    cluster_str <- paste("cluster", as.character(i - 1))
    output[cluster_str] <- color_vec[i]
  }
  output
}

extract_and_add_colors <- function(apotc_obj, plot_df) {
  color_hashmap <- pair_colors_to_hash(apotc_obj)
  color_vec <- plot_df[[1]] # 1 is label
  for (i in seq_along(color_vec)) {
    color_vec[i] <- color_hashmap[[color_vec[i]]]
  }
  plot_df %>% dplyr::mutate("color" = color_vec)
}

# in the future should probably make a fake ggplot cluster legend on the right
# side by inserting scatterplt points? (like in the seurat UMPA plot) But also,
# the problem is that it becomes inconsistent with the clone size legend :/

#' @noRd
#' @source https://stackoverflow.com/questions/649454/what-is-the-best-way-to-average-two-colors-that-define-a-linear-gradient#:~:text=NewColor%20%3D%20sqrt((R1%5E2%2BR2%5E2)/2)%2Csqrt((G1%5E2%2BG2%5E2)/2)%2Csqrt((B1%5E2%2BB2%5E2)/2)
get_average_hex <- function(hex1, hex2) {
  grDevices::rgb(
    t((col2rgb(hex1) + col2rgb(hex2)) / 2), maxColorValue = 255
  )
}

# for if hexes were square rooted first, typically from image files
get_average_image_hex_color <- function(hex1, hex2) {
  grDevices::rgb(
    t(sqrt(((col2rgb(hex1)^2) + (col2rgb(hex2)^2)) / 2)), maxColorValue = 255
  )
}
