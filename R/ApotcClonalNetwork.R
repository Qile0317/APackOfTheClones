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
#' @param clonesize_range integer vector of length 2. Sets the range of clone
#' sizes to keep when counting shared clones. The first element is the lower
#' bound (inclusive), and the second element is the upper bound (inclusive).
#' Defaults to `c(1L, Inf)`.
#' @param only_cluster integer vector indicating which clusters to keep when
#' counting shared clones. Note that this cannot conflict with
#' `exclude_cluster`.
#'
#' @return a named list where each name is a clonotype, each element is a
#' numeric indicating which seurat cluster(s) its in, in no particular order.
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

    clonesize_range = c(1L, Inf),
    only_cluster = NULL#,
    # only_between = NULL #TODO
    # TODO export format
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
        clone_size_lowerbound = clonesize_range[1],
        clone_size_upperbound = clonesize_range[2],
        included_cluster = create_cluster_truth_vector(
            only_cluster, exclude_cluster = NULL, get_num_clusters(apotc_obj)
        )
    )
}

getSharedClones_error_handler <- function() {
    args <- get_parent_func_args()
    check_apotc_identifiers(args)
    args$clonesize_range %>% typecheck(is_integer_pair)
    args$only_cluster %>% typecheck(is_an_integer, is_integer, is.null)
}

create_cluster_truth_vector <- function(
    only_cluster, exclude_cluster, num_clusters
) {
    if (is.null(only_cluster) && is.null(exclude_cluster)) {
        return(rep(TRUE, num_clusters))
    }

    if (!is.null(exclude_cluster)) {
        truth_indicies <- exclude_cluster
    } else {
        truth_indicies <- only_cluster
    }
    truth_val <- ifelse(!is.null(exclude_cluster), TRUE, FALSE)

    truth_vec <- rep(truth_val, num_clusters)
    truth_vec[as.integer(truth_indicies)] <- (!truth_val)
    truth_vec
}

# input: an ApotcData object
# output: a named list where each name is a clonotype, each element is a
# numeric indicating which seurat cluster(s) its in. If exclude_unique_clones,
# will filter out any clonotype with only length one.
# TODO - allow filtering based on original clone size
get_shared_clones <- function(
    apotc_obj,
    zero_indexed = FALSE,
    exclude_unique_clones = TRUE,
    clone_size_lowerbound = 1L,
    clone_size_upperbound = Inf,
    # TODO have heterogeneity bounds?
    included_cluster = TRUE
) {

    shared_clonotypes <- get_raw_shared_clones(apotc_obj)

    if (exclude_unique_clones) {
        shared_clonotypes <- remove_unique_clones(shared_clonotypes)
    }

    # TODO filter by clone size

    if (is_empty(shared_clonotypes) || all(included_cluster)) {
        return(shared_clonotypes)
    }

    filter_shared_clones_cluster(shared_clonotypes, included_cluster)
}

get_raw_shared_clones <- function(apotc_obj, zero_indexed = FALSE) {

    clonotype_map <- create_empty_int_hash(get_clonotypes(apotc_obj))
    clustered_clone_sizes <- get_raw_clone_sizes(apotc_obj, as_hash = TRUE)

    for (i in seq_along(clustered_clone_sizes)) {

        if (is_empty(clustered_clone_sizes[[i]])) next

        for (clonotype in names(clustered_clone_sizes[[i]])) {

            clonotype_map[[clonotype]] <- c(
                clonotype_map[[clonotype]], i - zero_indexed
            )

        }
    }

    # FIXME if some keys are too long, it becomes ... due to the environment implementation
    as.list(clonotype_map)
}

# FIXME
# - for single included cluster, should just connect from it to others
# - for multiple, only should linnk between them but also give option to expand out.
# - cirrently it expands out links - keep this implementation for if expand out is true.
filter_shared_clones_cluster <- function(shared_clonotypes, included_cluster) {
    results <- rcppFilterSharedClonesByClusterHelper(
        shared_clonotypes, names(shared_clonotypes), included_cluster
    )
    names(results[[2]]) <- results[[1]]
    results[[2]]
}

# overlay clone links on an APackOfTheClones plot
# maybe also make method to take in the plot directly?
# TODO - do some matrix visualization too, maybe use heatmap for clone sizes
overlay_shared_clone_links <- function(
    apotc_obj,
    result_plot,
    only_cluster,
    clonesize_range,
    link_type = "line", # TODO implement geom_ploygon link, also discuss way to make to better account for clonesize
    link_color_mode = "blend",
    link_alpha = 1,
    link_width = "auto",
    verbose = TRUE,
    link_mode = "default",
    extra_spacing = "auto" # not very relevant atm
) {
    shared_clones <- get_shared_clones(
        apotc_obj,
        zero_indexed = FALSE,
        exclude_unique_clones = TRUE,
        clone_size_lowerbound = clonesize_range[1],
        clone_size_upperbound = clonesize_range[2],
        create_cluster_truth_vector(
            only_cluster, exclude_cluster = NULL, get_num_clusters(apotc_obj)
        )
    )

    if (is_empty(shared_clones)) {
        if (verbose) message(
            "* no shared clonotypes with current filtering parameters"
        )
        return(result_plot)
    }
    
    if (identical(link_type, "line")) {
        link_dataframe <- compute_line_link_df( # FIXME
            apotc_obj, shared_clones, extra_spacing, link_mode
        )
    } else {
        stop(call. = FALSE, "no other link types are implemented yet")
    }

    link_dataframe <- add_link_colors(
        apotc_obj, link_dataframe, link_color_mode
    )

    # compute width and alpha (TODO make better formulas and figure out how to vectorize)
    link_alpha <- process_link_alpha(apotc_obj, link_alpha, verbose)
    link_width <- process_link_width(apotc_obj, link_width, verbose)

    result_plot %>%
        overlay_links(
            link_type = link_type, link_dataframe = link_dataframe,
            link_alpha = link_alpha, link_width = link_width
        )
}

# takes in a named list of clonotypes as names, the elements are numeric vectors
# indicating the seurat_cluster(s) they are in. If the numericvector is of length
# 1, remove the element. This is done in Rcpp to achieve true linear runtime.
remove_unique_clones <- function(shared_clonotypes) {
    results <- rcppRemoveUniqueClonesHelper(
        names(shared_clonotypes), shared_clonotypes
    )
    unique_clone_list <- results[[2]]
    names(unique_clone_list) <- results[[1]]
    unique_clone_list
}

# for link_mode = "line",
# outputs a dataframe with columns x1, x2, y1, y2, r1, r2
# the x and y correspond to line segments that originate and end at origins.
# the r corresponds to the radius of the circle.
# so will not be visually very appealing on its own.
compute_line_link_df <- function(
    apotc_obj, shared_clones, extra_spacing, link_mode
) {

    if (should_estimate(extra_spacing)) {
        extra_spacing <- 0 # TODO make better in future
    }

    if (link_mode != "default") {
        stop(call. = FALSE, "no other link modes are implemented yet")
    }

    rcppConstructLineLinkDf(
        clusterLists = get_clusterlists(apotc_obj),
        rawCloneSizes = get_raw_clone_sizes(apotc_obj),
        sharedClonotypeClusters = shared_clones,
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
                x = "x1", y = "y1", xend = "x2", yend = "y2",
                colour = "color"
            ),
            alpha = link_alpha,
            linewidth = link_width
        ) +
        ggplot2::scale_color_identity()
}
