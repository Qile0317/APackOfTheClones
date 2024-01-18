# input: apotc obj, clonal comparison function which defaults to equals.

# internal output1: 

# export-able output: a named list, each name is a clonotype, each value is which clonotype

# procedure:
# the easiest approach, taking advantage of sizes, since information is abstracted away -
# use a seed/any procedure to select clone sizes that match.

# intermediate flawed internal output (use Rcpp hashmap): list(clonotype=list(c1=s1,c2=s2,c3=s3))
# one alternative is each corresponding clone is linked exactly - each circle in the df_full_join should have the clone itself
# a simpler approach is to use any seed and link based on clonesize
# a final approach is to manually link and maybe re-adjust ordering of circles oneself.

# intermediate internal output dataframe of [x1 x2 y1 y2 r1 r2] - for future possibilities

# internal output for plotting LINES: [x1 x2 y1 y2] (subtract the radii and a rad decrease)
# - future: plot ribbons

# user API (S4 generic):
# should be an extra option in APOTCPlot `show_shared_clones`, `share_clone_line_spacing`, future:`shared_clone_linetype`


# internal master function for overlaying the clonal links based on mode: start with lines only.

overlay_shared_clone_links <- function(
    apotc_obj, result_plot, link_mode,
    link_type,
    extra_spacing, link_color_mode, link_alpha, link_width
) {
    shared_clones <- get_shared_clones(apotc_obj)

    if (link_type == "line") {
        link_dataframe <- compute_line_link_df(apotc_obj, shared_clones, link_mode)
    } else {
        stop(call. = FALSE, "no other link types are implemented yet")
    }

    link_dataframe

    # TODO overlay_links(hash::hash(as.list(environment())))
}

# input: an ApotcData object
# output: a named list where each name is a clonotype, each element is a
# numeric indicating which seurat cluster(s) its in. If exclude_unique_clones,
# will filter out any clonotype with only length one.
get_shared_clones <- function(apotc_obj, exclude_unique_clones = TRUE) {

    clonotype_map <- create_valueless_vector_hash(
        get_clonotypes(apotc_obj), numeric
    )

    num_unique_clonotypes <- length(clonotype_map)

    clustered_clonotypes <- lapply(get_raw_clone_sizes(apotc_obj), names)

    for (i in seq_along(clustered_clonotypes)) {
        if (is.null(clustered_clonotypes[[i]])) next
        for (clonotype in clustered_clonotypes[[i]]) {
            clonotype_map[[clonotype]] <- append(clonotype_map[[clonotype]], i)
        }
    }

    shared_clonotypes <- as.list(clonotype_map)
    if (!exclude_unique_clones) return(shared_clonotypes)
    remove_unique_clones(shared_clonotypes)
}

# takes in a named list of clonotypes as names, the elements are numeric vectors
# indicating the seurat_cluster(s) they are in. If the numericvector is of length
# 1, remove the element. This is done in Rcpp to achieve true linear runtime.
remove_unique_clones <- function(shared_clonotypes) {
    results <- rcppRemoveUniqueClonesHelper(names(shared_clonotypes), shared_clonotypes)
    unique_clone_list <- results[[2]]
    names(unique_clone_list) <- results[[1]]
    unique_clone_list
}

# TODO everything below is unfinished

compute_line_link_df <- function(apotc_obj, shared_clones, link_mode) {
    if (link_mode == "default") {
        return(rcppConstructLineLinkDf(
            clusterLists = get_clusterlists(apotc_obj),
            rawCloneSizes = get_raw_clone_sizes(apotc_obj),
            sharedClonotypeClusters = shared_clones
        ))
    } else {
        stop(call. = FALSE, "no other link modes are implemented yet")
    }
}

add_link_colors <- function(apotc_obj, link_dataframe, link_color_mode) {
    # TODO
    # - default - all the same color
    # - alt approach with lines: average color of two clones
}

# internal dispatch function to get a dataframe of line connections
# TODO should have exportable version with identifiers so user can get it and do their own thing

overlay_links <- function(args) {
    switch(args$link_type,
        "line" = return(overlay_line_links(args)),
        stop(call. = FALSE, "no other modes are implemented")
    )
}

overlay_line_links <- function(args) {
    args$result_plot + ggplot2::geom_segment(
        data = args$link_dataframe,
        mapping = ggplot2::aes(
            x = x1, y = y1, xend = x2, yend = y2 # TODO other params
        )
    )
}
