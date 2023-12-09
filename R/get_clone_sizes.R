# script probably will need to be revamped for multi-sample experiments. Might
# even end up involving heavy TCR/BCR sequence comparisons

# count the raw clone from the integrated seurat object, where col_condition is
# a dplyr condition to be evaluated. In the main function a string version
# can be mapped to conditions as shortcuts. assumes all inputs are valid
# returns the clone sizes list to be put in the ApotcData object
# TODO should make another function with a user wrapper AND as a getter
count_raw_clone_sizes <- function(
  seurat_obj, num_clusters, clonecall, metadata_filter_condition
) {
  if (!is.null(metadata_filter_condition)) {
    seurat_metadata <- seurat_obj@meta.data %>%
      dplyr::filter(substitute(metadata_filter_condition))
  } else {
    seurat_metadata <- seurat_obj@meta.data
  }

  clonotype_df <- data.frame(
    "clusters" = seurat_metadata[["seurat_clusters"]],
    "clonotype_ids" = seurat_metadata[[clonecall]]
  )

  # aggregate the raw counts
  freq_df <- stats::aggregate(
    as.formula("clonotype_ids ~ clusters"), clonotype_df, function(x) table(x)
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

# TODO convert input to metadata filter condition from string

get_processed_clone_sizes <- function(apotc_obj) {
  raw_tabled_clone_sizes <- apotc_obj@clone_sizes
  processed_sizes <- init_list(apotc_obj@num_clusters, list())
  for (i in 1:apotc_obj@num_clusters) {
    if (!is_empty_table(raw_tabled_clone_sizes[[i]][[1]])) {
      processed_sizes[[i]] <- apotc_obj@clone_scale_factor *
        sqrt(as.numeric(raw_tabled_clone_sizes[[i]][[1]]))
    }
  }
  processed_sizes
}
