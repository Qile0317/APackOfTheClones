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
  color_vec <- cluster_dataframe$label
  for (i in seq_along(color_vec)) {
    color_vec[i] <- color_hashmap[[color_vec[i]]]
  }
  cluster_dataframe %>% dplyr::mutate("color" = color_vec)
}

# new version that simply takes the readily existing colors in an seurat
# object, and adds the colors to the dataframe

# pair colors to hashmap
pair_colors_to_hash <- function(apotc_obj) {
  color_vec <- apotc_obj@cluster_colors
  output <- hash::hash()
  for (i in seq_len(get_num_clusters(apotc_obj))) {
    cluster_str <- paste("cluster", as.character(i - 1))
    output[cluster_str] <- color_vec[i]
  }
  output
}

extract_and_add_colors <- function(apotc_obj, plot_df) {
  color_hashmap <- pair_colors_to_hash(apotc_obj)
  color_vec <- plot_df$label
  for (i in seq_along(color_vec)) {
    color_vec[i] <- color_hashmap[[color_vec[i]]]
  }
  plot_df %>% dplyr::mutate("color" = color_vec)
}

#' @noRd
#' @source https://stackoverflow.com/questions/649454
get_average_hex <- function(hex1, hex2) {
  grDevices::rgb(
    t((grDevices::col2rgb(hex1) + grDevices::col2rgb(hex2)) / 2),
    maxColorValue = 255
  )
}

scale_hex_brightness <- function(hex, scaling_factor) {
  
  if (all(scaling_factor == 1)) return(hex)

  hsv_color <- hex %>%
    grDevices::col2rgb() %>%
    grDevices::rgb2hsv()

  hsv_color["v", ] <- bound_num(hsv_color["v", ] * scaling_factor, 0, 1)
  
  grDevices::hsv(hsv_color["h", ], hsv_color["s", ], hsv_color["v", ])
}
