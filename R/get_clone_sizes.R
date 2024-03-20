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
#' `scRepertoire`.
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
#' countCloneSizes(combined_pbmc, "nt", orig.ident = c("P17B", "P17L"))
#'
countCloneSizes <- function(
  seurat_obj, clonecall = "strict", extra_filter = NULL, ...
) {

    # type check inputs
    if (!is_seurat_object(seurat_obj))
        stop("`seurat_obj` must be a Seurat object.")
    typecheck(clonecall, is_a_character)
    typecheck(extra_filter, is.null, is_a_character)
    check_filtering_conditions(as.list(environment()), frame_level = 1)

    # setup variables
    clonecall <- .theCall(seurat_obj@meta.data, clonecall)
    filter_string <- parse_to_metadata_filter_str(
        metadata_filter = extra_filter, varargs_list = list(...)
    )

    if (is_valid_filter_str(filter_string)) { # TODO probably use seurat's built-in version :P, also probably shoudl allow for symbolic filtering.
        seurat_obj <- subsetSeuratMetaData(seurat_obj, filter_string)
    }

    count_raw_clone_sizes(
        seurat_obj = seurat_obj,
        num_clusters = get_num_total_clusters(seurat_obj),
        clonecall = clonecall
    )
}

# TODO create S4 generic to allow getting it from run_id, as an Apotc Getter

# count the raw clone from the integrated seurat object from the METADATA
# FIXME - for only 1 cluster doenst work
count_raw_clone_sizes <- function(
  seurat_obj, num_clusters, clonecall
) {

  # aggregate the raw counts
  freq_df <- stats::aggregate(
    stats::as.formula(paste(clonecall, "~ seurat_clusters")),
    seurat_obj@meta.data,
    table
  )

  # compile the tabled counts into a list of table objects
  clone_sizes <- init_list(num_clusters, create_empty_table())

  if (nrow(freq_df) == 0) return(clone_sizes)

  if (nrow(freq_df) == 1) {
    clone_sizes[[freq_df$seurat_clusters[1]]] <- as.table(freq_df[[2]][1, ])
    for (i in 1:num_clusters) names(dimnames(clone_sizes[[i]])) <- "" #?
    return(clone_sizes)
  }

  for (elem in enumerate(freq_df$seurat_clusters)) {
    clone_sizes[[val1(elem)]] <- freq_df[[2]][ind(elem)][[1]]
  }

  clone_sizes
}
