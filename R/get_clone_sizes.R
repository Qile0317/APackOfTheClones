# script probably will need to be revamped for multi-sample experiments. Might
# even end up involving heavy TCR/BCR sequence comparisons

# inefficient counting of raw clonotype sizes, and adding it to the apotc
add_raw_clone_sizes <- function(
    apotc_obj, integrated_seurat_obj, colname = "raw_clonotype_id"
) {
    df <- data.frame(
      "clusters" = integrated_seurat_obj@meta.data[["seurat_clusters"]],
      "clonotype_ids" = integrated_seurat_obj@meta.data[[colname]]
    )
    # aggregate the raw counts
    freq_df <- aggregate(clonotype_ids ~ clusters, df, function(x) table(x))

    # add the tabled counts, purposefully not modifying them
    cluster_indicies <- as.numeric(freq_df[[1]]) # converts to one based indexing!
    num_valid_clusters <- length(cluster_indicies)
    index <- 1
    for (i in 1:apotc_obj@num_clusters) {
        if (index > num_valid_clusters) {
            break
        }
        if (i != cluster_indicies[index]) {
            apotc_obj@clone_sizes[[i]] <- NA
        }
        apotc_obj@clone_sizes[[i]] <- freq_df[[2]][index]
        index <- index + 1
    }
    apotc_obj
}

get_transformed_clone_sizes <- function(
    sizelist, clone_scale_factor, num_clusters
) {
    output_sizes <- vector("list", num_clusters)
    for (i in 1:num_clusters) {
        curr_sizes <- sizelist[[i]]
        if (is.null(curr_sizes) || (!isnt_empty(curr_sizes))) {
            output_sizes[[i]] <- list()
            next
        }
        output_sizes[[i]] <- sqrt(curr_sizes) * clone_scale_factor
    }
    output_sizes
}

# might have to redo this
# multiply all tables in a raw clonotype freq list by a scale factor and sqrt (for the next update)
get_processed_clone_sizes <- function(apotc_obj) {
    raw_tabled_clone_sizes <- apotc_obj@clone_sizes
    processed_sizes <- vector("list", apotc_obj@num_clusters)
    for (i in 1:apotc_obj@num_clusters) {
      if (is.null(raw_tabled_clone_sizes[[i]][[1]])) {
        processed_sizes[[i]] <- list()
        next
      }
      processed_sizes[[i]] <- sqrt(as.numeric(raw_tabled_clone_sizes[[i]][[1]])) *
        apotc_obj@clone_scale_factor
    }
    processed_sizes
}

# memory and speed inefficient counting of clonotypes within each cluster,
# for the old version of the package
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
