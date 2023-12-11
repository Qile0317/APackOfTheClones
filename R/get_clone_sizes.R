# script probably will need to be revamped for multi-sample experiments. Might
# even end up involving heavy TCR/BCR sequence comparisons

# count the raw clone from the integrated seurat object from the METADATA
# TODO should make another function with a user wrapper AND as a getter
count_raw_clone_sizes <- function(
  seurat_obj, num_clusters, clonecall
) {

  # aggregate the raw counts
  freq_df <- stats::aggregate(
    as.formula(paste(clonecall, "~ seurat_clusters")),
    seurat_obj@meta.data,
    function(x) table(x)
  )

  # compile the tabled counts, purposefully not modifying them
  cluster_indicies <- as.numeric(freq_df[[1]]) # converts to one based indexing!
  num_valid_clusters <- length(cluster_indicies)
  index <- 1
  clone_sizes <- init_list(num_clusters, table(NULL))

  for (i in 1:num_clusters) {
    if (index > num_valid_clusters) {
      break
    }
    if (i != cluster_indicies[index]) {
      next
    }
    clone_sizes[[i]] <- freq_df[[2]][index]
    index <- index + 1
  }

  clone_sizes
}

get_processed_clone_sizes <- function(apotc_obj) {
  raw_tabled_clone_sizes <- apotc_obj@clone_sizes
  processed_sizes <- init_list(apotc_obj@num_clusters, list())
  for (i in 1:apotc_obj@num_clusters) {
    if (!is_empty_table(raw_tabled_clone_sizes[[i]])) {
      processed_sizes[[i]] <- apotc_obj@clone_scale_factor *
        sqrt(as.numeric(raw_tabled_clone_sizes[[i]][[1]]))
    }
  }
  processed_sizes
}
