# plotting for devs to test a single clusteer

plot_single_cluster_layout <- function(sizes) {
  a <- circle_layout(sizes)
  a <- df_full_join(list(a))
  plot_clusters(a)
}

plot_single_cluster <- function(cluster_list) {
  plot_clusters(df_full_join(list(cluster_list)))
}