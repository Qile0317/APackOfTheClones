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
#' clones so that only the clonotypes with overall `top` frequency/count is
#' returned. If the input is integer, it indicates to only include the top
#' clones overall by frequency ranking. If a numeric in (0, 1), this indicates
#' to include the `top` clones in that *proportion* of the frequencies.
#' @param top_per_cl integer or numeric in (0, 1) - if not null, filters the
#' output clones so that for each seurat cluster, only the clonotypes with the
#' `top_per_cl` frequency/count is preserved when aggregating shared clones.
#' If the input is integer, it indicates to only include the `top_per_cl` clones
#' overall by frequency ranking per cluster. If a numeric in (0, 1), this
#' indicates to include the `top_per_cl` clones in that *proportion* of the
#' frequencies for each cluster. Note that if inputted in conjunction with
#' `top`, it will get the *intersection* of the clonotypes filtered each way.
#'
#' @return a named list where each name is a clonotype, each element is a
#' numeric indicating which seurat cluster(s) its in, in no particular order.
#' If no shared clones are present, the output is an empty named list().
#' @export
#'
#' @examples
#' data("combined_pbmc")
#'
#' getSharedClones(combined_pbmc)
#'
#' getSharedClones(
#'     combined_pbmc,
#'     orig.ident = c("P17B", "P18B"), # this is a named subsetting parameter
#'     clonecall = "aa"
#' )
#'
#' getSharedClones(
#'     combined_pbmc,
#'     clonecall = "aa",
#'     top_per_cl = 0.5
#' )
#'
#' # doing a run and then getting the clones works too
#' combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "run1", verbose = FALSE)
#' getSharedClones(combined_pbmc, run_id = "run1")
#'
getSharedClones <- function(
    seurat_obj,

    reduction_base = "umap",
    clonecall = "strict",
    ...,
    extra_filter = NULL,
    run_id = NULL,

    top = NULL, # FIXME = 1 crashes session - probably due to no matches? test.
    top_per_cl = NULL

    # TODO top_shared and top_sh_per_cl. again intersect everything
) {
    # handle inputs
    varargs_list <- list(...)
	args <- hash::hash(as.list(environment()))
    getSharedClones_error_handler()

    # get the apotcdata
    apotc_obj <- getApotcDataIfExistsElseCreate(
        seurat_obj, infer_object_id_if_needed(args, varargs_list)
    )

    get_shared_clones(
        apotc_obj,
        zero_indexed = FALSE,
        exclude_unique_clones = TRUE,
        top_clones = top,
        top_per_cluster = top_per_cl
    )
}

getSharedClones_error_handler <- function() {
    args <- get_parent_func_args()
    check_apotc_identifiers(args)
    typecheck(args$top, is_an_integer, is_a_numeric_in_0_1, is.null)
    typecheck(args$top_per_cl, is_an_integer, is_a_numeric_in_0_1, is.null)
    # TODO
}

# input: an ApotcData object
# output: a named list where each name is a clonotype, each element is a
# numeric indicating which seurat cluster(s) its in. If exclude_unique_clones,
# will filter out any clonotype with only length one. (not shared)
# TODO - allow filtering based on original clone size
get_shared_clones <- function(
    apotc_obj,
    zero_indexed = FALSE,
    exclude_unique_clones = TRUE,
    top_per_cluster = NULL,
    top_clones = NULL #, # these are top in all so maybe theres no match - the top in shared clones is a diff story
    # clone_size_lowerbound = 1L,
    # clone_size_upperbound = Inf,
    # TODO have heterogeneity bounds?
) {
    apotc_obj %>%
        get_raw_clone_sizes() %>%
        filter_clonesizes_if_needed(top_clones, top_per_cluster) %>%
        get_raw_shared_clones(zero_indexed) %>%
        remove_unique_clones_if(exclude_unique_clones)
}

get_raw_shared_clones <- function(clustered_clone_sizes, zero_indexed = FALSE) {

    clonotype_map <- clustered_clone_sizes %>%
        unlist() %>%
        names() %>%
        create_empty_int_hash()

    for (i in seq_along(clustered_clone_sizes)) {

        if (is_empty(clustered_clone_sizes[[i]])) next

        for (clonotype in names(clustered_clone_sizes[[i]])) {

            clonotype_map[[clonotype]] <- c(
                clonotype_map[[clonotype]], i - zero_indexed
            )

        }
    }

    # FIXME if some keys are too long, it becomes ...
    as.list(clonotype_map)
}

filter_clonesizes_if_needed <- function(
    clone_sizes, top_clones, top_per_cluster
) {

    if (is.null(top_clones) && is.null(top_per_cluster)) {
        return(clone_sizes)
    }

    if (!is.null(top_clones)) {
        clone_sizes_by_top <- filter_top_clones(clone_sizes, top_clones)
    }

    if (!is.null(top_per_cluster)) {
        clone_sizes_by_top_per_cl <- filter_top_by_cluster(
            clone_sizes, top_per_cluster
        )
    }

    if (!is.null(top_clones)) {
        clone_sizes <- clone_sizes_by_top
        if (is.null(top_per_cluster)) return(clone_sizes)
    }

    intersect_common_table_lists(
        clone_sizes, clone_sizes_by_top_per_cl
    )
}

filter_top_clones <- function(clone_sizes, top_clones) {

    top_clonotypes <- get_top_clonotypes(clone_sizes, top_clones)

    lapply(clone_sizes, function(x) {

        if (is_empty(x)) return(x)

        filtered_positions <- which(names(x) %in% top_clonotypes)

        # this is for an edge case: (len(x) = 1) => (x[TRUE] != x) no idea why
        if ((length(filtered_positions) == 1) && filtered_positions == 1) {
            return(x)
        }

        filtered_x <- x[filtered_positions]
        if (is_empty(filtered_x)) create_empty_table() else filtered_x
    })
}

# will sort clones - maybe undesirable?
filter_top_by_cluster <- function(clone_sizes, top_clones) {

    clone_sizes <- sort_each_table(clone_sizes, desc = TRUE)

    if (is_a_numeric_in_0_1(top_clones)) {
        clone_table_handler <- function(x) {
            if (is_empty(x)) return(x)
            x[1:round(length(x) * top_clones)]
        }
    } else if (is_an_integer(top_clones)) {
        clone_table_handler <- function(x) {
            if (is_empty(x)) return(x)
            x[1:min(top_clones, length(x))]
        }
    }

    lapply(clone_sizes, clone_table_handler)

    # TODO
}

# overlay clone links on an APackOfTheClones plot
# maybe also make method to take in the plot directly?
# TODO - do some matrix visualization too, maybe use heatmap for clone sizes
overlay_shared_clone_links <- function(
    apotc_obj,
    shared_clones,
    result_plot,
    only_cluster, # if length 1, radiate from that cl. In future ver, between pairs
    link_type = "line", # TODO implement geom_ploygon link, also discuss way to make to better account for clonesize
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

# takes in a named list of clonotypes as names, the elements are numeric vectors
# indicating the seurat_cluster(s) they are in. If the numericvector is of length
# 1, remove the element. This is done in Rcpp to achieve true linear runtime.
remove_unique_clones_if <- function(shared_clonotypes, should_actually_remove) {

    if (!should_actually_remove) return(shared_clonotypes)

    results <- rcppRemoveUniqueClonesHelper(
        names(shared_clonotypes), shared_clonotypes
    )
    unique_clone_list <- results[[2]]
    names(unique_clone_list) <- results[[1]]
    unique_clone_list

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
        oneIndexedSourceClusterIndex = ifelse(is.null(only_cluster), -1, only_cluster),
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
# TODO should have exportable version with identifiers so user can get it and do their own thing

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
