#' @title
#' count the number of clonotype sizes in a seurat object combined with a
#' VDJ library overall or by cluster
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Get clonotype frequencies from a seurat object's metadata, either as one
#' whole table, or in a list of tables, based on the current / some custom
#' ident of each cell. Note that depending on the ident (indicated by the
#' `by_cluster` argument) there may be more or less clonotypes counted based
#' on the number of rows containing NA for that column of that ident if it
#' isn't the active ident.
#'
#' @param seurat_obj a seurat object combined with a VDJ library with
#' `scRepertoire`.
#' @inheritParams RunAPOTC
#' @param by_cluster Logical or Character. If `TRUE`, will output a list of
#' table objects, with the table at each index corresponding to level in
#' Idents(). Each table's names are the clonotype name indicated by `clonecall`
#' after filtering, while the values are the actual clone sizes. If `FALSE`,
#' outputs just the aggregate clone sizes for all cells. Note that if `FALSE`,
#' the output should be identical to that produced by
#' `mergeCloneSizes(countCloneSizes(..., by_cluster = TRUE))`. Otherwise, this
#' argument can also be a character indicating some column in the seurat object
#' metadata to use a cell identity guiding (e.g. `"seurat_clusters"`).
#' @param sort_decreasing a logical or NULL. If `TRUE`/`FALSE`, sorts each/the
#' table by clonotype frequency with largest/smallest clones first with a stable
#' sorting algorithm, and if NULL, no order is guaranteed but the output is
#' deterministic.
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

    seurat_obj <- set_meta_ident_col(
        seurat_obj, alt_ident = if_a_logical_convert_null(by_cluster)
    )

    if (is_valid_filter_str(filter_string)) {
        seurat_obj <- subsetSeuratMetaData(seurat_obj, filter_string)
    }

    seurat_obj <- ident_into_seurat_clusters(seurat_obj)

    clustered_clone_sizes <- count_raw_clone_sizes(
        seurat_obj = seurat_obj,
        ident_levels = get_ident_levels(seurat_obj, "seurat_clusters"),
        clonecall = clonecall,
        named = TRUE
    )

    if (is_false(by_cluster)) {
        return(mergeCloneSizes(clustered_clone_sizes, sort_decreasing))
    }

    if (!is.null(sort_decreasing)) {
        clustered_clone_sizes <- sort_each_clone_size_table(
            clustered_clone_sizes, sort_decreasing
        )
    }

    clustered_clone_sizes
}

countCloneSizes_arg_checker <- function() {
    args <- get_parent_func_args()

    if (!is_seurat_object(args$seurat_obj))
        stop("`seurat_obj` must be a Seurat object.")
    typecheck(args$clonecall, is_a_character)
    typecheck(args$extra_filter, is.null, is_a_character)
    check_filtering_conditions(as.list(environment()), frame_level = 2)
    typecheck(args$by_cluster, is_a_logical, is_a_character)
    typecheck(args$sort_decreasing, is_a_logical, is.null)
}

count_raw_clone_sizes <- function(
    seurat_obj, ident_levels, clonecall, named = FALSE
) {

    # get the na filtered clonotype metadata
    clonotype_cluster_df <- seurat_obj@meta.data %>%
        select_cols(clonecall, "seurat_clusters") %>%
        stats::na.omit()

    # initialize output clone sizes - named!
    clone_sizes <- init_empty_table_list(ident_levels)

    # if no valid rows, return empty clone sizes
    if(nrow(clonotype_cluster_df) == 0) return(clone_sizes %>% unname_if(!named))

    # aggregate the raw counts
    freq_df <- stats::aggregate(
        stats::as.formula(paste(clonecall, "~ seurat_clusters")),
        clonotype_cluster_df,
        table
    )

    if (nrow(freq_df) == 1) {

        if (is.matrix(freq_df[[2]])) {
            one_cluster_clone_table <- as_table(freq_df[[2]][1, ])
        } else { # edge case with only one cluster and one clonotype
            one_cluster_clone_table <- freq_df[[2]]
            names(one_cluster_clone_table) <- clonotype_cluster_df[1, clonecall]
            one_cluster_clone_table <- as_table(one_cluster_clone_table)
        }

        clone_sizes[[
            as.character(freq_df$seurat_clusters[1])
        ]] <- one_cluster_clone_table

        return(clone_sizes %>% unname_if(!named))
    }

    for (ident in as.character(freq_df[[1]])) {
        clone_sizes[[ident]] <- freq_df[[2]][
            as.character(freq_df[[1]]) == ident][[1]]
    }

    clone_sizes %>% unname_if(!named)
}

#' @title
#' Merge a list of Clustered Clonotype Frequency Tables
#'
#' @description
#' The list of clustered clonotype frequencies from [countCloneSizes]
#' can be merged by this function to a frequency table of all clonotypes
#' similar to the data that can be seen in the seurat object metadata.
#' By default, this function sorts the table with largest clonotypes first,
#' and this may be useful for quickly gauging which clonotypes are the most
#' expanded overall.
#'
#' @param clustered_clone_sizes the output of [countCloneSizes].
#' @param sort_decreasing a logical or NULL. If `TRUE`/`FALSE`, sorts the
#' table by clonotype frequency with largest/smallest clones first, and if
#' NULL, no order is guaranteed but the output is deterministic.
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

    if (is_list_of_empty_tables(clustered_clone_sizes)) {
        warning("no clones are present")
        return(clustered_clone_sizes)
    }

    clustered_clone_sizes %>%
        aggregate_clone_sizes(sort_decreasing) %>%
        convert_named_numeric_to_table()

}

# union the output of count_raw_clone_sizes to a named numeric
# (not table but exactly like one)
# depending on the ident, there may be more or less NA rows dropped when aggregating by cluster
# the following outputs 317
# (countCloneSizes(combined_pbmc, by_cluster = "mito.genes") %>% mergeCloneSizes() == (countCloneSizes(combined_pbmc,by_cluster=FALSE))) %>% sum
# if "mito.genes" is replaced with TRUE, output is 341
#
aggregate_clone_sizes <- function(
    clone_sizes, sort_decreasing = NULL, top_clones = NULL
) {

    if (!is.null(top_clones)) sort_decreasing <- TRUE # technically not needed
    union_clone_sizes <- union_list_of_tables(clone_sizes, sort_decreasing)

    if (is.null(top_clones) || is_empty(union_clone_sizes)) {
        return(union_clone_sizes)
    }
   
    filter_top_clonesize(union_clone_sizes, top_clones)
}

# helper for filtering aggregate_clone_sizes
filter_top_clonesize <- function(union_clone_sizes, top_clones) {

    unique_clone_sizes_desc <- union_clone_sizes %>%
        get_unique_table_values()

    num_unique_sizes <- length(unique_clone_sizes_desc)

    clonesize_lowerbound <- unique_clone_sizes_desc[
        ifelse(
            test = is_a_numeric_in_0_1(top_clones),
            yes = round(num_unique_sizes * top_clones),
            no = min(top_clones, num_unique_sizes)
        )
    ]
   
    union_clone_sizes[union_clone_sizes >= clonesize_lowerbound]
}

get_top_clonotypes <- function(clone_sizes, top_clones) {
    names(aggregate_clone_sizes(clone_sizes, top_clones = top_clones))
}

sort_each_clone_size_table <- function(x, decreasing) {
    sort_each_table(x, decreasing)
}

## abstract higher order functions for filtering clustered clone sizes

# filter clone sizes with two 2 arg functions, the first function
# takes in 2 args, clone_sizes, and filter_1_arg, and so on for the
# second one. If either arg is null, the clone size will not be
# filtered fwith the corresponding function. If both are non-null,
# the results are set intersected
filter_clonesize_2way_if_need <- function(
    clone_sizes,
    filter_func_1, filter1_arg,
    filter_func_2, filter2_arg
) {
    if (is.null(filter1_arg) && is.null(filter2_arg)) {
        return(clone_sizes)
    }

    if (!is.null(filter1_arg) && is.null(filter2_arg)) {
        return(filter_func_1(clone_sizes, filter1_arg))
    }

    if (is.null(filter1_arg) && !is.null(filter2_arg)) {
        return(filter_func_2(clone_sizes, filter2_arg))
    }

    intersect_common_table_lists(
        filter_func_1(clone_sizes, filter1_arg),
        filter_func_2(clone_sizes, filter2_arg)
    )
}
