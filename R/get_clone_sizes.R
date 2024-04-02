#' @title
#' count the number of clonotype sizes in a seurat object combined with a
#' VDJ library overall or by cluster
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Get clonotype frequencies from a seurat object's meta.data slot.
#'
#' @param seurat_obj a seurat object combined with a VDJ library with the
#' `scRepertoire`.
#' @inheritParams RunAPOTC
#' @param by_cluster If `TRUE`, will output a list of table objects, with the
#' table at each index corresponding to each seurat cluster index. Each table's
#' names are the clonotype name indicated by `clonecall` after filtering, while
#' the values are the actual clone sizes. Else, outputs just the aggregate clone
#' sizes for all cells. Note that if `FALSE`, the output should be identical to
#' that produced by `mergeCloneSizes(countCloneSizes(..., by_cluster = TRUE))`
#' @param sort_decreasing a logical or NULL. If `TRUE`/`FALSE`, sorts the
#' table by clonotype frequency with largest/smallest clones first, and if
#' NULL, no order is guaranteed but the output is deterministic.
#'
#' @return A list of tables or a single table depending on `by_cluster`
#' @export
#' 
#' @seealso [mergeCloneSizes]
#'
#' @examples
#' data("combined_pbmc")
#'
#' countCloneSizes(combined_pbmc)
#' countCloneSizes(combined_pbmc, "aa")
#' countCloneSizes(combined_pbmc, "nt", orig.ident = c("P17B", "P17L"))
#'
countCloneSizes <- function(
    seurat_obj,
    clonecall = "strict",
    extra_filter = NULL,
    ...,
    by_cluster = TRUE,
    sort_decreasing = NULL
) {
    countCloneSizes_arg_checker()

    # setup variables
    clonecall <- .theCall(seurat_obj@meta.data, clonecall)
    filter_string <- parse_to_metadata_filter_str(
        metadata_filter = extra_filter, varargs_list = list(...)
    )

    if (is_valid_filter_str(filter_string)) { # TODO probably use seurat's built-in version :P, also probably shoudl allow for symbolic filtering.
        seurat_obj <- subsetSeuratMetaData(seurat_obj, filter_string)
    }

    clustered_clone_sizes <- count_raw_clone_sizes(
        seurat_obj = seurat_obj,
        num_clusters = get_num_total_clusters(seurat_obj),
        clonecall = clonecall
    )

    if (!by_cluster) {
        return(mergeCloneSizes(clustered_clone_sizes, sort_decreasing))
    }

    if (!is.null(sort_decreasing)) {
        clustered_clone_sizes <- lapply(
            clustered_clone_sizes,
            function(x) sort(x, decreasing = sort_decreasing, method = "radix")
        )
    }

    clustered_clone_sizes
}

# TODO create S4 generic to allow getting it from run_id, as an Apotc Getter

countCloneSizes_arg_checker <- function() {
    args <- get_parent_func_args()

    if (!is_seurat_object(args$seurat_obj))
        stop("`seurat_obj` must be a Seurat object.")
    typecheck(args$clonecall, is_a_character)
    typecheck(args$extra_filter, is.null, is_a_character)
    check_filtering_conditions(as.list(environment()), frame_level = 2)
    typecheck(args$by_cluster, is_a_logical)
    typecheck(args$sort_decreasing, is_a_logical, is.null)
}

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

#' @title
#' Merge a list of Clustered Clonotype Frequency Tables
#'
#' @description
#' The list of clustered clonotype frequencies from [countCloneSizes]
#' can be merged by this function to a frequency table of all clonotypes
#' similar to the data that can be seen in the seurat object metadata.
#' By default, this function sorts the table with largest clonotypes first,
#' and this may be useful for quickly guaging which clonotypes are the most
#' expanded overall.
#'
#' @param clustered_clone_sizes the output of [countCloneSizes].
#' @inheritParams countCloneSizes
#'
#' @return a table object
#' @export
#'
#' @seealso [countCloneSizes]
#'
#' @examples
#' clustered_clone_sizes <- countCloneSizes(get(data("combined_pbmc")))
#' mergeCloneSizes(clustered_clone_sizes)
#'
mergeCloneSizes <- function(clustered_clone_sizes, sort_decreasing = TRUE) {

    typecheck(clustered_clone_sizes, is_output_of_countCloneSizes)
    typecheck(sort_decreasing, is_a_logical, is.null)

    if (is_empty(clustered_clone_sizes)) {
        warning("no clone are present")
        return(clustered_clone_sizes)
    }

    clustered_clone_sizes %>%
        aggregate_clone_sizes(sort_decreasing) %>%
        convert_named_numeric_to_table()

}

# union the output of count_raw_clone_sizes to a named numeric
# (not table but exactly like one)
aggregate_clone_sizes <- function(
    clone_sizes, sort_decreasing = NULL, top_clones = NULL
) {

    if (!is.null(top_clones)) sort_decreasing <- TRUE
    union_clone_sizes <- union_list_of_tables(clone_sizes, sort_decreasing)

    if (is.null(top_clones)) return(union_clone_sizes)
    num_clones <- length(union_clone_sizes)

    if (is_an_integer(top_clones)) {
        return(union_clone_sizes[1:min(top_clones, num_clones)])
    }

    if (is_a_numeric_in_0_1(top_clones)) {
        return(union_clone_sizes[1:round(num_clones * top_clones)])
    }

    # TODO more filtering

}

get_top_clonotypes <- function(clone_sizes, top_clones) {
    names(aggregate_clone_sizes(clone_sizes, top_clones = top_clones))
}
