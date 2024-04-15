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
#' should be based on.
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
getReductionCentroids <- function(seurat_obj, reduction) {
  get_cluster_centroids(
    seurat_obj = seurat_obj,
    reduction = user_get_reduc_obj(seurat_obj, reduction),
    get_ident_levels(seurat_obj)
  )
}

user_get_reduc_obj <- function(seurat_obj, reduction) {
    if (!is_seurat_object(seurat_obj))
        stop(call. = FALSE, "`seurat_obj` not a seurat object!")
    if (!is_a_character(reduction))
        stop(call. = FALSE, "`reduction` must be one character")
    seurat_obj@reductions[[attempt_correction(seurat_obj, reduction)]]
}

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

# data initialization utils

init_numeric <- function(x, init_val = 0) {

    if (is_an_integer(x)) {
        named_vals <- numeric(x)
    } else {
        named_vals <- structure(numeric(length(x)), names = x)
    }

    sapply(named_vals, function(y) init_val)
}

# table/(list of table) utils

create_empty_table <- function() {
    structure(
        integer(0),
        dim = 0L,
        dimnames = structure(list(NULL), names = ""),
        class = "table"
    )
}

convert_table_to_named_numeric <- function(x) {
    structure(as.numeric(x), names = names(x))
}

as_table <- function(x) ensure_proper_table(as.table(x))

ensure_proper_table <- function(x) {
    if (is.null(names(dimnames(x)))) names(dimnames(x)) <- ""
    x
}

strip_to_numeric <- function(table_obj) {
    convert_table_to_named_numeric(table_obj)
}

convert_named_numeric_to_table <- function(named_numeric_vector) {
    output <- as.integer(named_numeric_vector)
    names(output) <- names(named_numeric_vector)
    output <- as.table(output)
    names(dimnames(output)) <- ""
    output
}

# intersect tables, ***ASSUMING*** that for elements with common names,
# they have the exact same value (frequency).
intersect_common_tables <- function(t1, t2) {
    t1[names(t1) %in% intersect(names(t1), names(t2))]
}

# vectorized version of intersect_common_tables, assuming lists are the
# same length, and ignores if either table is empty
intersect_common_table_lists <- function(l1, l2) {
    operate_on_same_length_lists(intersect_common_tables, l1, l2)
}

union_list_of_tables <- function(x, sort_decreasing = NULL, as_table = FALSE) {

    if (is_list_of_empty_tables(x)) {
        if (as_table) return(create_empty_table())
        return(numeric(0))
    }

    if (length(x) == 1) {
        return(sort_if_desc_arg_not_null(
            convert_table_to_named_numeric(x[[1]]), sort_decreasing
        ))
    }

    x <- x %>%
        lapply(convert_table_to_named_numeric) %>%
        union_list_of_named_numerics()

	if (!is.null(sort_decreasing)) {
        x <- sort(x, decreasing = sort_decreasing, method = "radix")
    }
	if (as_table) x <- convert_named_numeric_to_table(x)
    x
}

union_list_of_named_numerics <- function(x) {
    rcppUnionListOfNamedNumericsHelper(x)
}

sort_each_table <- function(x, desc = FALSE) {
    lapply(x, function(y) {
        if (is_empty(y)) return(y)
        as_table(sort(y, decreasing = desc, method = "radix"))
    })
}

init_empty_table_list <- function(x) {
    init_list(x, create_empty_table())
} 

is_list_of_empty_tables <- function(x) is_empty(x) || all(sapply(x, is_empty))

# dataframe utils

select_cols <- function(df, ...) df[unlist(list(...))]

# hash::hash utilities

create_hash_from_keys <- function(keys, init_vals = NULL) {
    keys <- unique(keys)
    numkeys <- length(keys)
    switch(as.character(numkeys),
        "0" = hash::hash(),
        "1" = hash::hash(keys, init_vals),
        hash::hash(keys, init_list(numkeys, init_vals))
    )
}

create_valueless_vector_hash <- function(key, vector_type) {
    create_hash_from_keys(key, vector_type(0))
}

create_empty_int_hash <- function(keys) {
    create_valueless_vector_hash(keys, integer)
}

hash_from_tablelist <- function(tablelist) {
    lapply(
        tablelist,
        function(x) {
            if (!is_empty_table(x)) return(hash::hash(x))
            hash::hash()
        }
    )
}

# readability functions

sort_if_desc_arg_not_null <- function(x, desc) {
    if (is.null(desc)) return(x)
    sort(x, decreasing = desc, method = "radix")
}

sort_if_needed <- function(x, do_sort, desc = FALSE) {
    if (is.null(do_sort) || is_false(do_sort)) return(x)
    sort(x, decreasing = desc)
}

if_a_logical_convert_null <- function(x) {
    if (!is_a_logical(x)) return(x)
    return(NULL)
}

contains_duplicates <- function(v) anyDuplicated(v) != 0

is_false <- function(x) identical(x, FALSE)

is_empty <- function(inp) (length(inp) == 0L) || identical(inp, hash::hash())

isnt_empty <- function(inp) !is_empty(inp)

isnt_na <- function(inp) !any(is.na(inp))

isnt_empty_nor_na <- function(inp) isnt_empty(inp) && isnt_na(inp)

is_empty_table <- is_empty

is_int <- function(num) all(num == as.integer(num))

should_estimate <- function(obj, auto_str = "auto") identical(obj, auto_str)

should_assume <- should_estimate

should_change <- function(obj) !is.null(obj)

should_compute <- function(x) is.null(x)

as_expression <- function(...) {
    parse(text = paste0(unlist(list(...)), collapse = ""))
}

subset_dataframe <- function(df, filter_string) {
    df %>% dplyr::filter(eval(as_expression(filter_string)))
}

# ggplot2 utils

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

name_latest_layer <- function(plt, new_name) {
    if (is.null(names(plt$layers))) {
        names(plt$layers) <- rep("", length(plt$layers))
    }
    names(plt$layers)[length(plt$layers)] <- new_name
    plt
}

remove_ggplot_layers <- function(ggplot_obj, layer_indices) {
  ggplot_obj$layers[layer_indices] <- NULL
  ggplot_obj
}

get_ggplot_data <- function(x) x$data

set_ggplot_data <- function(ggplot_obj, new_data) {
  ggplot_obj$data <- new_data
  ggplot_obj
}

# naming utils

secretly_init_name <- function(x) {
    names(x) <- rep("", length(x))
    x
}

unname_if_empty <- function(l) if (is_empty(l)) unname(l) else l

unname_if <- function(x, do_unname) if (do_unname) unname(x) else x

# math utils

bound_num <- function(num, lowerbound, upperbound) {
    min(max(num, lowerbound), upperbound)
}

is_bound_between <- function(num, lowerbound, upperbound) {
    (num >= lowerbound) & (num <= upperbound)
}

add <- function(x, y) x + y
subtract <- function(x, y) x - y
multiply <- function(x, y) x * y
divide <- function(x, y) x / y

is_even <- function(x) x %% 2 == 0
is_odd <- function(x) x %% 2 == 1

# special mutation operators from StackOverflow

create_mutator <- function(binary_operator) {
    function(var, val) {
        eval(
            call("<-", substitute(var), binary_operator(var, val)),
            envir = parent.frame()
        )
    }
}

"%+=%" <- create_mutator(add)
"%*=%" <- create_mutator(multiply)

# iteration utils

get_unique_pairs_up_to <- function(x) {
    if (x <= 1) return(list())
    rcppGetUniquePairsUpTo(as.integer(x), oneIndexed = TRUE)
}

zip <- function(...) mapply(list, ..., SIMPLIFY = FALSE)

enumerate <- function(..., zero_indexed = FALSE) {
    zip(seq_along(..1) - zero_indexed, ...)
}
ind <- function(elem) elem[[1]]
val <- function(elem, index) elem[[index + 1]]
val1 <- function(elem) val(elem, 1)

# spelling utils

is_vowel <- function(s) tolower(s) %in% c("a", "e", "i", "o", "u")
starts_with_vowel <- function(s) is_vowel(get_first_char(s))

prepend_indefinite_article <- function(s, exclude = c("NULL")) {
    if (s %in% exclude) return(s)
    paste("a", ifelse(starts_with_vowel(s), "n", ""), " ", s, sep = "")
}

get_right_of_dollarsign <- function(x) sapply(strsplit(x, "\\$"), getlast)

get_first_char <- function(s) substr(s, 1, 1)

strip_spaces <- function(s) gsub(" ", "", s)

strip_and_lower <- function(s) strip_spaces(tolower(s))

strip_unquoted_spaces <- function(input_str) {

    all_parts <- strsplit(input_str, "'")

    for (i in seq_along(input_str)) {
    parts <- all_parts[[i]]

        for (j in seq_along(parts)) {
            if (is_odd(j)) parts[j] <- strip_spaces(parts[j])
        }

        input_str[i] <- Reduce(function(...) paste(..., sep = "'"), parts)

        if (is_even(length(parts))) {
            input_str[i] <- paste(input_str[i], "'", sep = "")
        }
    }

    input_str
}

user_attempt_correction <- function(
    s,
    strset,
    stop_msg_start,
    modifiers = list(tolower, trimws, strip_unquoted_spaces, strip_spaces)
) {

    # check if the string is already present in strset and if yes return
    match_indices <- which(s == strset)
    if (length(match_indices) == 1) return(s)

    get_only_similar_word_or_null <- function(modifier) {
        match_indices <- which(modifier(s) == modifier(strset))
        if (length(match_indices) != 1) return(NULL)
        message(paste(
            "* assuming `", s, "` corresponds to `",
            strset[match_indices], "`", sep = ""
        ))
        strset[match_indices]
    }

    for (modifier in append(identity, modifiers)) {
        potential_unique_similar_word <- get_only_similar_word_or_null(modifier)
        if (!is.null(potential_unique_similar_word)) {
            return(potential_unique_similar_word)
        }
    }

    for (ij in get_unique_pairs_up_to(length(modifiers))) {
        potential_unique_similar_word <- get_only_similar_word_or_null(
            modifier = function(x) modifiers[[ij[1]]](modifiers[[ij[2]]](x))
        )
        if (!is.null(potential_unique_similar_word)) {
            return(potential_unique_similar_word)
        }
    }
    
    stop(
        stop_msg_start, " `", s, "`, did you mean: `",
        closest_word(s, strset), "`?",
        call. = FALSE
    )
}

closest_word <- function(s, strset) {
    strset <- unique(strset)
    if (length(strset) == 1) return(strset)

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

init_list <- function(x, init_val = NULL) {
    l <- vector("list", ifelse(is_an_integer(x), x, length(x)))
    if (!is_an_integer(x)) names(l) <- as.character(x)
    if (is.null(init_val)) return(l)
    lapply(l, function(x) init_val)
}

getlast <- function(x, n = 1) {
    index <- length(x) - n + 1
    if (is.list(x)) return(x[[index]])
    x[index]
}

# get first non empty element in a list, assuming it exists
get_first_nonempty <- function(lst) {
    for (el in lst) if (isnt_empty(el)) return(el)
}

# operate on non-empty elements of two lists of the same length
# with a 2-argument function
operate_on_same_length_lists <- function(func, l1, l2) {
    l <- init_list(length(l1), list())
    for (i in seq_along(l1)) {
        if (is_empty(l1[[i]]) || is_empty(l2[[i]])) next
        l[[i]] <- func(l1[[i]], l2[[i]])
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

repr_as_string <- function(input) {
    insertchar <- ""
    if (is.character(input)) insertchar <- "'"
    to_string_rep_with_insert(v = input, insert = insertchar)
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

# Seurat utils

subsetSeuratMetaData <- function(
    seurat_obj, filter_string, error_param = "extra_filter"
) {
	seurat_obj@meta.data <- subset_dataframe(seurat_obj@meta.data, filter_string)

	if (nrow(seurat_obj@meta.data) == 0) {
		stop(call. = FALSE, paste(
			"please check `", error_param, "`, ",
			"no rows in the seurat metadata match the filter condition",
            sep = ""
		))
	}

	seurat_obj
}

count_clones <- function(seurat_obj, clonecall) {
  sum(!is.na(seurat_obj@meta.data[[clonecall]]))
}

# TODO check if identical with ApotcData order
get_ident_levels <- function(seurat_obj, custom_ident = NULL) {

    if (is.null(custom_ident)) {
        return(as.character(levels(seurat_obj@active.ident)))
    }

    if (!is.null(seurat_obj@meta.data[["__active.ident__"]])) {
        return(as.character(levels(seurat_obj@meta.data[["__active.ident__"]])))
    }

    as.character(levels(as.factor(seurat_obj@meta.data[[custom_ident]])))
}

get_num_total_clusters <- function(seurat_obj) {
  length(levels(seurat_obj@meta.data[["seurat_clusters"]]))
}

# seurat reduction related functions

any_reduction_exists <- function(seurat_obj) {
    reduction_names <- get_curr_reduc_names(seurat_obj)
    !(is.null(reduction_names) || identical(reduction_names, character(0)))
}

get_curr_reduc_names <- function(seurat_obj) {
    names(seurat_obj@reductions)
}

get_2d_embedding <- function(seurat_obj, reduction) {
    SeuratObject::Embeddings(object = seurat_obj, reduction = reduction)[, 1:2]
}

attempt_correction <- function(seurat_obj, reduction) {

    if (!any_reduction_exists(seurat_obj)) {
        stop("No dimensional reductions detected")
    }

    curr_reductions <- get_curr_reduc_names(seurat_obj)

    if (identical(strip_and_lower(reduction), "t-sne")) {
        if (!any("t-sne" == strip_and_lower(curr_reductions))) {
            reduction <- "tsne"
        }
    }

    user_attempt_correction(
      reduction,
      strset = curr_reductions,
      stop_msg_start = "Invalid reduction"
    )
}
