# This dummy function definition is included with the package to ensure that
# 'tools::package_native_routine_registration_skeleton()' generates the required
# registration info for the 'run_testthat_tests' symbol
(function() {
    .Call("run_testthat_tests", FALSE, PACKAGE = "APackOfTheClones")
})

# progress bar functions

progress_bar <- function (x = 0, max = 100) {
    percent <- 100 * (x / max)
    cat(sprintf(
        '\r[%-50s] %d%%',
        paste(rep('=', percent * 0.5), collapse = ''),
        floor(percent)
    ))
}

start_progress_bar <- function(verbose = TRUE) {
    if (verbose) {
        progress_bar(0, 1)
    }
}

end_progress_bar <- function(verbose = TRUE) {
    if (verbose) {
        progress_bar(1, 1)
    }
}

print_completion_time <- function(start_time, digits = 3, newline = FALSE) {
    end_time <- Sys.time()
    if (newline) cat("\n")
    message(paste(
        "\nCompleted successfully, time elapsed:",
        round(as.numeric(end_time - start_time), digits),
        "seconds\n"
    ))
}

# readability functions

create_empty_table <- function() {
    structure(
        integer(0),
        dim = 0L,
        dimnames = structure(list(NULL), names = ""),
        class = "table"
    )
}

is_empty <- function(inp) identical(inp, list())
isnt_empty <- function(inp) !identical(inp, list())

isnt_na <- function(inp) !any(is.na(inp))

isnt_empty_nor_na <- function(inp) isnt_empty(inp) && isnt_na(inp)

is_empty_table <- function(inp) identical(inp, table(NULL))

is_int <- function(num) all(num == as.integer(num))

should_estimate <- function(obj, auto_str = "auto") identical(obj, auto_str)

should_assume <- should_estimate

should_change <- function(obj) !is.null(obj)

should_compute <- function(x) is.null(x)

strip_spaces <- function(s) gsub(" ", "", s)

strip_and_lower <- function(s) strip_spaces(tolower(s))

strip_unquoted_spaces <- function(input_str) {
  parts <- strsplit(input_str, "'")[[1]]
  for (i in seq_along(parts)) {
    if (is_odd(i)) parts[i] <- strip_spaces(parts[[i]])
  }
  
  joined_str <- Reduce(function(...) paste(..., sep = "'"), parts)

  if (is_even(length(parts))) {
    return(paste(joined_str, "'", sep = ""))
  }
  joined_str
}

#' @title Get the xmin, xmax, ymin, ymax of a ggplot object
#' @return list(xr = c(xmin, xmax), yr = c(ymin, ymax))
#' @noRd
get_plot_dims <- function(plt) {
    built_plt_layout <- ggplot2::ggplot_build(plt)$layout
    list(
        xr = built_plt_layout$panel_scales_x[[1]]$range$range,
        yr = built_plt_layout$panel_scales_y[[1]]$range$range
    )
}

get_xr <- function(p) {
    if (ggplot2::is.ggplot(p)) {
        return(ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$range$range)
    }
    p[[1]]
}

get_yr <- function(p) {
    if (ggplot2::is.ggplot(p)) {
        return(ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    }
    p[[2]]
}

is_seurat_object <- function(obj) inherits(obj, "Seurat")

is_sce_object <- function(obj) inherits(obj, "SingleCellExperiment")

is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")

is_seurat_or_sce_object <- function(obj) {
    is_seurat_object(obj) || is_sce_object(obj)
}

is_a_character <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.character(x)
}

is_an_integer <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.integer(x)
}

is_a_numeric <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.numeric(x)
}

is_a_logical <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.logical(x)
}

# math utils

bound_num <- function(num, lowerbound, upperbound) {
    min(max(num, lowerbound), upperbound)
}

is_bound_between <- function(num, lowerbound, upperbound) {
    num >= lowerbound && num <= upperbound
}

add <- function(x, y) x + y
subtract <- function(x, y) x - y

is_even <- function(x) x %% 2 == 0
is_odd <- function(x) x %% 2 == 1

# spelling related functions

user_attempt_correction <- function(s, strset, stop_msg_start) {
    if (any(strip_and_lower(s) == strip_and_lower(strset))) {
        return(s)
    }
    stop(paste(
        stop_msg_start,
        ": '", s, "', did you mean: ",
        closest_word(s, strset),
        sep = ""
    ))
}

closest_word <- function(s, strset = c("umap", "tsne", "pca")) {
    strset <- unique(strset)
    strset_lowercase <- tolower(strset)
    s <- tolower(s)

    closest_w <- strset_lowercase[1]
    closest_dist <- utils::adist(s, closest_w)
    for(i in 2:length(strset_lowercase)) {
        curr_dist <- utils::adist(s, strset_lowercase[i])
        if (curr_dist < closest_dist) {
            closest_w <- strset[i]
            closest_dist <- curr_dist
        }
    }
    closest_w
}

# list utilities

init_list <- function(num_elements, init_val = NULL) {
    l <- vector("list", num_elements)
    for (i in 1:num_elements) {
        l[[i]] <- init_val
    }
    l
}

getlast <- function(x) UseMethod("getlast")
getlast.default <- function(x) x[length(x)]
getlast.list <- function(x) x[[length(x)]]

# operate on non-empty elements of two lists of the same length
# with a 2-argument function
operate_on_same_length_lists <- function(func, l1, l2) {
    l <- init_list(length(l1), list())
    for (i in seq_along(l1)) {
        if (isnt_empty(l1[[i]]) && isnt_empty(l2[[i]])) {
            if (!(is.null(l1[i]) || is.null(l2[i]))) {
                l[[i]] <- func(l1[[i]], l2[[i]])
            }
        }
    }
    l
}

move_coord_list_by_same_amount <- function(
    coord_list, original_coord_list, new_coord_list
) {
    operate_on_same_length_lists(
        func = add,
        l1 = coord_list,
        l2 = operate_on_same_length_lists(
            func = subtract,
            l1 = new_coord_list,
            l2 = original_coord_list
        )
    )
}

#' Take a list of character vectors and join each element of the vectors
#' together, separating each character by sep. Currently recursive which
#' will be bad for larger inputs :P
#' @return a character vector
#' @noRd
construct_prefix_vector <- function(params, sep = "_") {
    unlist(join_list_of_characters(params, sep))
}

join_list_of_characters <- function(params, sep = "_") {

    if (length(params) == 2) {
        l2 <- params[[2]]
    } else {
        l2 <- construct_prefix_vector(params[2:length(params)])
    }

    operate_on_same_length_lists(
        func = function(x, y) paste(x, y, sep = sep),
        l1 = params[[1]],
        l2 = l2
    )
}

# S3 method to represent vectors as strings

repr_as_string <- function(input, ...) {
    UseMethod("repr_as_string")
}

repr_as_string.character <- function(input, ...) {
    to_string_rep_with_insert(v = input, insert = "'")
}

repr_as_string.default <- function(input, ...) {
    to_string_rep_with_insert(v = input, insert = "")
}

# represent vector as string - doesnt take into account of names!
to_string_rep_with_insert <- function(v, insert) {
    if (length(v) == 1) {
        return(paste(insert, v, insert, sep = ""))
    }

    output <- ""
    for (x in v) {
        output <- paste(output, insert, x, insert, ",", sep = "")
    }
    paste("c(", substr(output, 1, nchar(output) - 1), ")", sep = "")
}

repr_as_string.list <- function(input, ...) {
    list_names <- names(input)
    named <- !is.null(list_names)
    output <- ""

    for (i in seq_along(input)) {
        if (named && list_names[i] != "")
            output <- paste(output, list_names[i], "=", sep = "")
        output <- paste(output, repr_as_string(input[[i]]), ",", sep = "")
    }
    paste("list(", substr(output, 1, nchar(output) - 1), ")", sep = "")
}

# Seurat utils

#' @title
#' Calculate seurat cluster centroids based on a Dimensional reduction
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Utility function to calculate the physical xy coordinates of each seurat
#' cluster based on a dimensional reduction already present in the object.
#' The results are returned in a list with the length of the number of distinct
#' seurat clusters based on the seurat_obj `meta.data`.
#'
#' @param seurat_obj input seurat object with the dimensional reduction of
#' choice already present, and seurat clusters computed.
#' @param reduction character. The reduction that the centroid calculation
#' should be based on. Currently, can only be "umap", "tsne", or "pca".
#'
#' @return A list of the length of the number of distinct clusters in the
#' seurat object metadata, where each element of the list is a numeric vector
#' of length 2, with the numbers corresponding to the x and y coordinate
#' respectively of the seurat cluster with the corresponding index.
#'
#' @export
#'
#' @examples
#' data("combined_pbmc")
#' getReductionCentroids(combined_pbmc, reduction = "umap")
#'
getReductionCentroids <- function(seurat_obj, reduction = "umap") {
  get_cluster_centroids(
    seurat_obj = seurat_obj,
    reduction = user_get_reduc_obj(seurat_obj, reduction),
    passed_in_reduc_obj = TRUE
  )
}

# Returns the number of valid barcodes that are not NA's
count_tcr_barcodes <- function(seurat_obj) {
  sum(!is.na(seurat_obj@meta.data[["barcode"]]))
}

count_clones <- function(seurat_obj, clonecall) {
  sum(!is.na(seurat_obj@meta.data[[clonecall]]))
}

# get the percent of NA's in the metadata barcode column for the message
percent_na <- function(seurat_obj) {
  num_barcodes <- length(seurat_obj@meta.data[["barcode"]])
  100 * (num_barcodes - count_tcr_barcodes(seurat_obj)) / num_barcodes
}

get_rna_assay_barcodes <- function(seurat_obj) {
    seurat_obj@assays[["RNA"]]@data@Dimnames[[2]]
}

# seurat cluster related functions

count_num_clusters <- function(seurat_obj) {
  data.table::uniqueN((seurat_obj@meta.data[["seurat_clusters"]]))
}

get_num_total_clusters <- function(seurat_obj) {
  length(levels(seurat_obj@meta.data[["seurat_clusters"]]))
}

# seurat reduction related functions

get_2d_embedding <- function(seurat_obj, reduction) {
  seurat_obj@reductions[[reduction]]@cell.embeddings[, 1:2]
}

attempt_correction <- function(s) {

    s <- strip_and_lower(s)
    s <- ifelse(test = identical(s, "t-sne"), yes = "tsne", no = s)

    user_attempt_correction(
      s,
      strset = c("umap", "tsne", "pca"),
      stop_msg_start = "Invalid reduction"
    )
}

is_reduction_name <- function(reduction) {
  any(reduction == c("umap", "tsne", "pca"))
}

user_get_reduc_obj <- function(seurat_obj, reduction) {
  if (!is_seurat_object(seurat_obj)) {
    stop("not a seurat object!")
  }

  reduction <- attempt_correction(reduction)

  if (!is_reduction_name(reduction)) {
    stop(paste("A", reduction, "has not been ran on the seurat object"))
  }

  seurat_obj@reductions[[reduction]]
}
