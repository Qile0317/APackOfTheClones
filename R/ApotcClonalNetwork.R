# TODO matrix/heatmap visualization for number of shared clones
# TODO somehow take into account sizes of shared clones

#' @title Compute a list of clonotypes that are shared between seurat clusters
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is a convenience function for users to get a list of clonotypes
#' which are shared between seurat clusters. If `run_id` is inputted, then the
#' function will attempt to get the shared clonotypes from the corresponding
#' APackOfTheClones run generated from [RunAPOTC]. Otherwise, it will use the
#' filtering / subsetting parameters to generate the shared clones by
#' internally running [RunAPOTC] with the same parameters and returning the
#' shared clones.
#'
#' @inheritParams RunAPOTC
#' @param top integer or numeric in (0, 1) - if not null, filters the output
#' clones so that only the shared clonotypes with counts the top `top` count /
#' proportion (for numeric in (0, 1) input) shared clones are kept.
#' @param top_per_cl integer or numeric in (0, 1) - if not null, filters the
#' output clones so that for each seurat cluster, only the clonotypes with the
#' `top_per_cl` frequency/count is preserved when aggregating shared clones,
#' in the same way as the above. Note that if inputted in conjunction with
#' `top`, it will get the *intersection* of the clonotypes filtered each way.
#' @param intop integer or numeric in (0, 1) - if not null, filters the raw
#' clone sizes before computing the shared clonotypes so that only the
#' clonotypes that have their overall size in the top `intop` largest sizes
#' (if it is integer, else the `intop` proportion) are kept. To emphasize,
#' this argument ***does not necessarily*** return the `top` shared clones
#' and likely a little less, because this filters the raw clone sizes, of
#' which, its very likely that not all those clones end up being shared.
#' @param intop_per_cl integer or numeric in (0, 1) - if not null, filters
#' the raw *clustered* clone sizes before computing shared clones, so that
#' for every clone in a seurat cluster, the top `intop_per_cl` count /
#' proportion (for numeric in (0, 1) input) clones are kept.
#'
#' @return a named list where each name is a clonotype, each element is a
#' numeric indicating which seurat cluster(s) its in, in no particular order.
#' If no shared clones are present, the output is an empty list.
#' @export
#'
#' @examples
#' data("combined_pbmc")
#'
#' getSharedClones(combined_pbmc)
#'
#' getSharedClones(
#'     combined_pbmc,
#'     orig.ident = c("P17B", "P18B"), # a named subsetting parameter
#'     clonecall = "aa"
#' )
#'
#' # extract shared clones from a past RunAPOTC run
#' combined_pbmc <- RunAPOTC(
#'     combined_pbmc, run_id = "foo", verbose = FALSE
#' )
#'
#' getSharedClones(
#'     combined_pbmc, run_id = "foo", top = 5
#' )
#'
#' # doing a run and then getting the clones works too
#' combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "run1", verbose = FALSE)
#' getSharedClones(combined_pbmc, run_id = "run1")
#'
getSharedClones <- function( # TODO FIXME output now can not longer just be numbers - but convert to numbers if the original factor levels are 1:num_clusters
    seurat_obj,
    reduction_base = "umap",
    clonecall = "strict",
    ...,
    extra_filter = NULL,
    alt_ident = NULL,
    run_id = NULL,

    top = NULL,
    top_per_cl = NULL,

    intop = NULL,
    intop_per_cl = NULL
) {
    # handle inputs
    varargs_list <- list(...)
    getSharedClones_error_handler()

    get_shared_clones(
        apotc_obj = getApotcDataIfExistsElseCreate(
            seurat_obj, run_id, environment(), ...
        ),
        zero_indexed = FALSE,
        exclude_unique_clones = TRUE,
        in_top_clones = intop,
        in_top_per_cluster = intop_per_cl,
        top_shared = top,
        top_shared_per_cluster = top_per_cl
    )
}

getSharedClones_error_handler <- function() {
    args <- get_parent_func_args()
    check_apotc_identifiers(args)
    typecheck(args$intop,
        is_a_positive_integer, is_a_numeric_in_0_1, is.null)
    typecheck(args$intop_per_cl,
        is_a_positive_integer, is_a_numeric_in_0_1, is.null)
    typecheck(args$top,
        is_a_positive_integer, is_a_numeric_in_0_1, is.null)
    typecheck(args$top_per_cl,
        is_a_positive_integer, is_a_numeric_in_0_1, is.null)
}

# input: an ApotcData object
# output: a named list where each name is a clonotype, each element is a
# numeric indicating which seurat cluster(s) its in. If exclude_unique_clones,
# will filter out any clonotype with only length one. (not shared)
# TODO FIXME make it based on ident
get_shared_clones <- function(
    apotc_obj,
    zero_indexed = FALSE,
    exclude_unique_clones = TRUE,
    in_top_clones = NULL,
    in_top_per_cluster = NULL,
    top_shared = NULL,
    top_shared_per_cluster = NULL#,
    # min_publicity = 0,
    # min_size = NULL,
    # min_size_per_cluster = NULL
) {

    raw_clone_sizes <- get_raw_clone_sizes(apotc_obj)

    raw_clone_sizes %>%
        filter_top_sizes_if_needed(in_top_clones, in_top_per_cluster) %>%
        # filter_min_clones_if_needed(min_size, min_size_per_cluster) %>% # again need to union
        get_raw_shared_clones(zero_indexed) %>%
        remove_unique_clones_if(exclude_unique_clones) %>%
        filter_shared_if_needed(
            raw_clone_sizes, top_shared, top_shared_per_cluster
        ) %>%
        unname_if_empty()
}

filter_top_sizes_if_needed <- function(
    clone_sizes, top_clones, top_per_cluster
) {
    clone_sizes %>% filter_clonesize_2way_if_need(
        filter_top_clones, top_clones,
        filter_top_by_cluster, top_per_cluster
    )
}

filter_top_clones <- function(clone_sizes, top_clones) {
    top_clonotypes <- get_top_clonotypes(clone_sizes, top_clones)
    clone_sizes %>% lapply(function(x) {
        if (is_empty(x)) x else as_table(x[names(x) %in% top_clonotypes])
    })
}

# will sort clones - could make it not sort but shouldnt matter
filter_top_by_cluster <- function(clone_sizes, top_clones) {

    if (is_list_of_empty_tables(clone_sizes)) return(clone_sizes)

    clone_sizes <- sort_each_table(clone_sizes, desc = TRUE)

    if (is_an_integer(top_clones)) {
        clone_table_handler <- function(x) {
            as_table(x[1:min(top_clones, length(x))])
        }
    }
    else if (is_a_numeric_in_0_1(top_clones)) {
        clone_table_handler <- function(x) {
            top_index <- round(length(x) * top_clones)
            if (top_index == 0) return(create_empty_table())
            as_table(x[1:top_index])
        }
    }

    lapply(
        clone_sizes, function(x) if (is_empty(x)) x else clone_table_handler(x)
    )
}

# # TODO
# filter_min_clones_if_needed <- function(
#     clone_sizes, min_size, min_size_per_cluster
# ) {
#     clone_sizes %>% filter_clonesize_2way_if_need(

#     )
# }

# filter_min_clones <- function(clone_sizes, min_size) { # integer or float in (0,1)

# }

get_raw_shared_clones <- function(clustered_clone_sizes, zero_indexed = FALSE) {

    if (is_list_of_empty_tables(clustered_clone_sizes)) return(list())

    shared_clonotype_map <- clustered_clone_sizes %>%
        unlist() %>%
        names() %>%
        create_empty_int_hash()

    for (i in seq_along(clustered_clone_sizes)) {

        if (is_empty(clustered_clone_sizes[[i]])) next

        for (clonotype in names(clustered_clone_sizes[[i]])) {

            shared_clonotype_map[[clonotype]] <- c(
                shared_clonotype_map[[clonotype]], i - zero_indexed
            )

        }
    }

    # FIXME if some keys are too long, it becomes ...
    as.list(shared_clonotype_map)
}

# takes in a named list of clonotypes as names, the elements are numeric vectors
# indicating the seurat_cluster(s) they are in. If the numericvector is of length
# 1, remove the element. This is done in Rcpp to achieve true linear runtime.
remove_unique_clones_if <- function(shared_clonotypes, should_actually_remove) {

    if (!should_actually_remove || is_empty(shared_clonotypes)) {
        return(shared_clonotypes)
    }

    results <- rcppRemoveUniqueClonesHelper(
        names(shared_clonotypes), shared_clonotypes
    )
    unique_clone_list <- results[[2]]
    names(unique_clone_list) <- results[[1]]
    unique_clone_list

}

filter_shared_if_needed <- function(
    shared_clonotypes, clone_sizes, top, top_per_cl
) {
    shared_clonotypes %>%
        filter_top_shared(clone_sizes, top) %>%
        filter_top_shared_per_cl(clone_sizes, top_per_cl)
}

filter_top_shared <- function(shared_clones, raw_clone_sizes, top) {
    
    if (is.null(top) || is_empty(shared_clones)) return(shared_clones)

    filter_top_shared_w_clone_table(
        shared_clones,
        mergeCloneSizes(raw_clone_sizes, sort_decreasing = NULL),
        top
    )

}

filter_top_shared_w_clone_table <- function(shared_clones, clone_table, top) {

    if (is_list_of_empty_tables(clone_table) || is_empty(shared_clones)) {
        return(shared_clones)
    }

    sorted_shared_clones <- sort(
        clone_table[names(shared_clones)], decreasing = TRUE, method = "radix"
    )

    if (is_an_integer(top)) {
        top_index <- min(length(shared_clones), top)
    } else if (is_a_numeric_in_0_1(top)) {
        top_index <- length(shared_clones) * top
    } else {
        stop("dev error: wrong variable type")
    }

    filtered_clones <- names(sorted_shared_clones)[1:top_index]
    filter_shared_that_contain(shared_clones, filtered_clones)
}

filter_top_shared_per_cl <- function(shared_clones, clone_sizes, top) {

    if (is.null(top) ||
        is_list_of_empty_tables(clone_sizes) ||
        is_empty(shared_clones)
    ) {
        return(shared_clones)
    }

    filtered_clones <- sort_each_table(clone_sizes, desc = TRUE) %>%
        lapply(function(x) {
            if (is_empty(x)) return(list())
            names(filter_top_shared_w_clone_table(shared_clones, x, top))
        }) %>%
        unlist()

    filter_shared_that_contain(shared_clones, filtered_clones)
}

filter_shared_that_contain <- function(shared_clones, filter_clonotypes) {
    shared_clones[names(shared_clones) %in% filter_clonotypes]
}

# overlay clone links on an APackOfTheClones plot
# maybe also make method to take in the plot directly?
overlay_shared_clone_links <- function(
    apotc_obj,
    shared_clones,
    result_plot,
    only_cluster, # TODO allow between pairs
    link_type = "line",
    link_color_mode = "blend",
    link_alpha = 1,
    link_width = "auto",
    verbose = TRUE,
    link_mode = "default",
    extra_spacing = "auto" # not very relevant atm
) {

    if (identical(link_type, "line")) {
        link_dataframe <- compute_line_link_df(
            apotc_obj, shared_clones, extra_spacing, link_mode, only_cluster
        )
    } else {
        stop(call. = FALSE, "no other link types are implemented yet")
    }

    link_dataframe <- add_link_colors(
        apotc_obj, link_dataframe, link_color_mode
    )

    result_plot %>%
        overlay_links(
            link_type = link_type,
            link_dataframe = link_dataframe,
            link_alpha = process_link_alpha(apotc_obj, link_alpha, verbose),
            link_width = process_link_width(apotc_obj, link_width, verbose)
        )
}

compute_line_link_df <- function(
    apotc_obj, shared_clones, extra_spacing, link_mode, only_cluster
) {

    if (link_mode != "default") {
        stop(call. = FALSE, "dev error: no other link modes are implemented")
    }

    if (should_estimate(extra_spacing)) {
        extra_spacing <- 0 # TODO make better in future
    }

    rcppConstructLineLinkDf(
        clusterLists = get_clusterlists(apotc_obj),
        rawCloneSizes = get_raw_clone_sizes(apotc_obj),
        sharedClonotypeClusters = shared_clones,
        oneIndexedSourceClusterIndex = ifelse(
            is.null(only_cluster), -1, only_cluster
        ),
        extraSpacing = extra_spacing - get_rad_decrease(apotc_obj)
    )
}

add_link_colors <- function(apotc_obj, link_dataframe, link_color_mode) {
    switch(link_color_mode,
        "blend" = return(add_blend_link_colors(apotc_obj, link_dataframe)),
        return(add_plain_link_colors(link_dataframe, link_color_mode))
    )
}

add_blend_link_colors <- function(apotc_obj, link_dataframe) {
    colors <- get_cluster_colors(apotc_obj)
    # extremeley cursed hack fix to pass R CMD check:
    eval(as_expression(
        "link_dataframe %>% dplyr::mutate(",
            "color = get_average_hex(colors[c1], colors[c2])",
        ")"
    ))
}

add_plain_link_colors <- function(link_dataframe, link_color) {
    link_dataframe$color <- rep(link_color, nrow(link_dataframe))
    link_dataframe
}

process_link_alpha <- function(apotc_obj, link_alpha, verbose) {
    if (should_estimate(link_alpha)) {
        link_alpha <- estimate_link_alpha(apotc_obj)
        if (verbose) message(
            "* setting `clone_link_alpha` to ", round(link_alpha, digits = 3)
        )
    }
    link_alpha
}

estimate_link_alpha <- function(apotc_obj) {
    0.75 # TODO make better in the future
}

process_link_width <- function(apotc_obj, link_width, verbose) {
    if (should_estimate(link_width)) {
        link_width <- estimate_link_width(apotc_obj)
        if (verbose) message(
            "* setting `clone_link_width` to ", round(link_width, digits = 3)
        )
    }
    link_width
}

estimate_link_width <- function(apotc_obj) {
    3 * get_clone_scale_factor(apotc_obj) #TODO improve
}

# internal dispatch function to get a dataframe of line connections
overlay_links <- function(
    result_plot, link_type, link_dataframe, link_alpha, link_width
) {
    switch(link_type,
        "line" = overlay_line_links(
            result_plot,
            link_dataframe,
            link_alpha,
            link_width
        ) %>% return()
        # should not get to any other case, this is just here for future extensions
    )
}

overlay_line_links <- function(
    result_plot, link_dataframe, link_alpha, link_width
) {
    result_plot +
        ggplot2::geom_segment(
            data = link_dataframe,
            mapping = apotc_aes_string(
                x = "x1",
                y = "y1",
                xend = "x2",
                yend = "y2",
                colour = "color"
            ),
            alpha = link_alpha,
            linewidth = link_width
        ) +
        ggplot2::scale_color_identity()
}
