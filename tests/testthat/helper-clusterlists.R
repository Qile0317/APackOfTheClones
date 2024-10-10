expect_apotc_data_clusterlist_equals_expected <- function(
    apotc_obj, expected_clusterlist, tolerance = 1e-6
) {
    testthat::expect_true(
        areGeometricallyEqualClusterLists(
            a = get_clusterlists(apotc_obj),
            b = expected_clusterlist,
            rad_decrease = get_rad_decrease(apotc_obj),
            clone_scale_factor = get_clone_scale_factor(apotc_obj),
            tolerance = tolerance
        )
    )
}

applyListAsArgsTo <- function(arglist, f, ...) do.call(f, arglist, ...)

areGeometricallyEqualClusterLists <- function(
    a, b, rad_decrease, clone_scale_factor, tolerance = 1e-6
) {

    if (get_num_clones(a) != get_num_clones(b)) {
        return(FALSE)
    }

    if (!setequal(get_radii(a), get_radii(b))) {
        return(FALSE)
    }

    if (!setequal(get_clonotypes(a), get_clonotypes(b))) {
        return(FALSE)
    }

    list_of_normed_ab <- list(a, b) %>%
        lapply(function(x)
            normalizeClusterListRadii(x, rad_decrease, clone_scale_factor)
        )
    
    list_of_normed_ab %>%
        append(list(threshold = tolerance)) %>%
        applyListAsArgsTo(areGeometricallyEqualNormalizedClusterLists)

}

normalizeClusterListRadii <- function(
    clusterlist, rad_decrease, clone_scale_factor
) {
    set_radii(
        clusterlist,
        (get_radii(clusterlist) + rad_decrease) / clone_scale_factor
    )
}

# compare two clusterlists that have their radii transformed back into their
# original clone counts.
#
# this assumes that the number of clones are equal, the radii are setequal,
# and that the overall clonotypes are setequal.
#
# returns TRUE if for all radii groups, the clonotypes are setequal
# (ignoring x, y) and the (x, y) pairs are setequal. If true, it means
# that the circle packing cluster made are geometrically identical and
# represent the same clones.
areGeometricallyEqualNormalizedClusterLists <- function(a, b, threshold = 1e-6) {

    list_of_clusterlist_dfs <- list(a, b) %>% lapply(function(clusterlist) {
        clusterlist %>%
            convert_to_dataframe("placeholder") %>%
            dplyr::mutate(r = r / min(r)) %>%
            dplyr::select(-label)
    })

    for (radius in unique(get_radii(a))) {
        
        clusterlist_dfs_filtered_by_curr_rad <- list_of_clusterlist_dfs %>%
            lapply(function(x) dplyr::filter(x, r = rad))

        are_dim_equal <- clusterlist_dfs_filtered_by_curr_rad %>%
            applyListAsArgsTo(function(a, b) identical(dim(a), dim(b)))
        if (!are_dim_equal) return(FALSE)
        
        clonotypes_are_setequal <- clusterlist_dfs_filtered_by_curr_rad %>%
            lapply(function(x) x["clonotype"]) %>%
            applyListAsArgsTo(setequal)

        if (!clonotypes_are_setequal) return(FALSE)

        xy_are_setequal <- clusterlist_dfs_filtered_by_curr_rad %>%
            lapply(function(a) dplyr::arrange(a, x, y)) %>%
            applyListAsArgsTo(function(df1, df2) {
                sum(df1 - df2) < ncol(df1) * nrow(df1) * threshold
            })

        if (!xy_are_setequal) return(FALSE)
    }

    return(TRUE)

}
