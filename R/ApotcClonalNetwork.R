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
    run_id = NULL
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
        exclude_unique_clones = TRUE
    )
}

getSharedClones_error_handler <- function() {
    args <- get_parent_func_args()
    check_apotc_identifiers(args)
}

# for the future - for two bit arrays, get A union neg B.
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
# will filter out any clonotype with only length one. (not shared)
# TODO - allow filtering based on original clone size
get_shared_clones <- function(
    apotc_obj,
    zero_indexed = FALSE,
    exclude_unique_clones = TRUE,
    top_per_cluster = NULL,
    top_clones = NULL # top clones = 1 creates some segfault that crashes R: not compatible with requested type list, int
    # clone_size_lowerbound = 1L,
    # clone_size_upperbound = Inf,
    # TODO have heterogeneity bounds?
) {

    clone_sizes <- get_raw_clone_sizes(apotc_obj)

    if (!is.null(top_clones)) {
        clone_sizes <- filter_top_clones(clone_sizes, top_clones)
    }

    if (!is.null(top_per_cluster)) {
        clone_sizes <- filter_top_by_cluster(clone_sizes, top_per_cluster)
    }

    shared_clonotypes <- get_raw_shared_clones(clone_sizes, zero_indexed)

    if (exclude_unique_clones) {
        shared_clonotypes <- remove_unique_clones(shared_clonotypes)
    }

    shared_clonotypes
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

# TODO FIXME segfaults for top_clones=1
filter_top_clones <- function(clone_sizes, top_clones) {

    top_clonotypes <- get_top_clonotypes(clone_sizes, top_clones)

    lapply(clone_sizes, function(x) {
        if (is_empty(x)) return(x)
        filtered_x <- x[names(x) %in% top_clonotypes]
        if (is_empty(filtered_x)) {
            return(create_empty_table())
        }
        filtered_x
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
    shared_clones <- get_shared_clones(
        apotc_obj,
        zero_indexed = FALSE,
        exclude_unique_clones = TRUE
    )

    if (is_empty(shared_clones)) {
        if (verbose) message(
            "* no shared clonotypes with current filtering parameters"
        )
        return(result_plot)
    }
    
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
remove_unique_clones <- function(shared_clonotypes) {
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
