#' @title count the number of clonotype sizes per cell cluster in a seurat
#' object combined with a VDJ library
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Get clonotype frequencies from a seurat object's meta.data slot.
#'
#' @inheritParams RunAPOTC
#' @param seurat_obj a seurat object combined with a VDJ library with the
#' `scRepertoire` workflow and [combineSeuratExpression].
#'
#' @return A list of table objects, with the table at each index corresponding
#' to each cluster index. Each table's names are the clonotype name indicated
#' by `clonecall` after filtering, while the values are the actual clone sizes.
#' @export
#'
#' @examples
#' data("combined_pbmc")
#'
#' countCloneSizes(combined_pbmc)
#' countCloneSizes(combined_pbmc, "aa")
#' countCloneSizes(combined_pbmc, "nts", orig.ident = c("P17B", "P17L"))
#'
countCloneSizes <- function(
  seurat_obj, clonecall = "strict", extra_filter = NULL, ...
) {

  # check inputs
  if (!is_seurat_object(seurat_obj))
    stop("`seurat_obj` must be a Seurat object.")
  if (!is_a_character(clonecall))
		stop("`clonecall` must be a character of length 1.")
  if (!is.null(extra_filter) && !is_a_character(extra_filter))
		stop("`extra_filter` must be a character of length 1.")
  check_filtering_conditions(as.list(environment()))

  clonecall <- .theCall(seurat_obj@meta.data, clonecall)
  filter_string <- parse_to_metadata_filter_str(
    metadata_filter = extra_filter, varargs_list = list(...)
  )

  if (is_valid_filter_str(filter_string)) {
    seurat_obj <- subsetSeuratMetaData(seurat_obj, filter_string)
  }

  count_raw_clone_sizes(
    seurat_obj = seurat_obj,
    num_clusters = get_num_total_clusters(seurat_obj),
    clonecall = clonecall
  )
}

# count the raw clone from the integrated seurat object from the METADATA
# TODO should make another function with a user wrapper AND as a getter
count_raw_clone_sizes <- function(
  seurat_obj, num_clusters, clonecall
) {

  # aggregate the raw counts
  freq_df <- stats::aggregate(
    stats::as.formula(paste(clonecall, "~ seurat_clusters")),
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
