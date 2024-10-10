#script for functions to deal with centroids and cluster coords

# All clusters of circles are "clusterlist" ADTs:
# [["x"]] numeric vector of x coordinates of all circles
# [["y"]] numeric vector of y coordinates of all circles
# [["rad"]] numeric vector of radii of all circles
# [["centroid"]] numeric vector of the cluster centroid x and y coordinate
# [["clRad"]] approximate radius of the cluster
# [["clonotype"]] clonotype based on original clonecall

# getters for a single clusterlist

get_x_coords <- function(l) l[[1]]
get_y_coords <- function(l) l[[2]]
get_radii <- function(l) l[[3]]
get_centroid <- function(l) l[[4]]
get_centroid_x <- function(l) get_centroid(l)[1]
get_centroid_y <- function(l) get_centroid(l)[2]
get_cluster_radius <- function(l) l[[5]]
get_clonotypes <- function(l) l$clonotype # not index due to legacy clusterlists

get_num_clones <- function(l) length(get_x_coords(l))

contains_clonotypes <- function(x) !is.null(get_clonotypes(x))

# setters for a single clusterlist
set_x_coords <- function(l, v) {l$x <- v; l}
set_y_coords <- function(l, v) {l$y <- v; l}
set_radii <- function(l, v) {l$rad <- v; l}
set_centroid <- function(l, v) {l$centroid <- v; l}
set_cluster_radius <- function(l, v) {l$clRad <- v; l}
set_clonotypes <- function(l, v) {l$clonotype <- v; l}

# convert clusterlist to dataframe, assuming its ***valid***
# returns a dataframe with columns:
# - label
# - x
# - y
# - r
# - clonotype
convert_to_dataframe <- function(
    clstr_list, seurat_cluster_index, detail = TRUE
) {
    if (!detail) {
      return(data.frame(
        "label" = paste("cluster", seurat_cluster_index),
        "x" = get_centroid_x(clstr_list),
        "y" = get_centroid_y(clstr_list),
        "r" = get_cluster_radius(clstr_list)
      ))
    }

    data.frame(
        "label" = rep(
            paste("cluster", seurat_cluster_index),
            get_num_clones(clstr_list)
        ),
        "x" = get_x_coords(clstr_list),
        "y" = get_y_coords(clstr_list),
        "r" = get_radii(clstr_list)
    ) %>%
        update_clusterlist_df(clstr_list)
}

update_clusterlist_df <- function(df, clusterlist) {

    if (contains_clonotypes(df) || !contains_clonotypes(clusterlist)) return(df)
   
    if (nrow(df) == 0) {
        df$clonotype <- character(0)
        return(df)
    }

    df %>% dplyr::mutate(clonotype = get_clonotypes(clusterlist))

}

# centroid finder for a matrix of [x, y, cluster]
find_centroids <- function(xyc_df, ident_levels) {

  cll <- split(xyc_df, factor(xyc_df[, 3])) %>%
    lapply(function(x) c(mean(x[, 1]), mean(x[, 2])))
   
  list_output <- init_list(ident_levels, list())

  for (ident_level in ident_levels) {
    if (is.null(cll[[ident_level]])) next
    list_output[[ident_level]] <- cll[[ident_level]]
  }

  unname(list_output)
}

# get reduction centroids from seurat obj, where the barcodes in the reduction
# cell embeddings will be filtered to be exactly the same as those left in the
# metadata incase it was additionally filtered.
get_cluster_centroids <- function(seurat_obj, reduction, ident_levels) {

  if (!is_a_character(reduction)) {
    reduc_coords <- reduction@cell.embeddings[, 1:2]
  } else {
    reduc_coords <- get_2d_embedding(seurat_obj, reduction)
  }

  find_centroids(
    xyc_df = data.frame(
      rcppFilterReductionCoords(
        seuratBarcodes = rownames(seurat_obj@meta.data),
        reductionCoords = reduc_coords
      ),
      "ident" = seurat_obj@meta.data[["seurat_clusters"]]
    ),
    ident_levels
  )
}

# if new coord is null, assumes current coords are from the
# cluster being centered at (0, 0), and transforms x and y by
# centroid. Otherwise, translates by new_coord.
trans_coord <- function(cluster, new_coord = NULL) {

  if (!is.null(new_coord)) {
    trans_vec <- new_coord
    cluster <- set_centroid(cluster, get_centroid(cluster) + new_coord)
  } else {
    trans_vec <- get_centroid(cluster)
  }

  cluster %>%
    set_x_coords(get_x_coords(cluster) + trans_vec[1]) %>%
    set_y_coords(get_y_coords(cluster) + trans_vec[2])

}

# MOVE clusterlist to a new centroid, irrespective of previous centroid
move_cluster <- function(cluster, new_coord) {
  dx <- cluster[[4]][1] - new_coord[1]
  dy <- cluster[[4]][2] - new_coord[2]

  cluster[[1]] <- cluster[[1]] - dx
  cluster[[2]] <- cluster[[2]] - dy
  cluster[[4]] <- new_coord
  cluster
}

# function to GET a list of centroids from a list of clusterlists,
read_centroids <- function(list_of_clusterlists) {
  lapply(list_of_clusterlists, function(x) {
    if (is_empty(x)) list() else get_centroid(x)
  })
}

areGeometricallyEqualClusterLists <- function(
    a, b, rad_decrease, clone_scale_factor, tolerance = 1e-6, verbose = FALSE
) {

    logIfVerbose <- function(...) if (verbose) message(glue(...))
    logSetDiffIfVerbose <- function(...) {
        logIfVerbose(getSetDiffAsListInStr(...))
    }

    if (get_num_clones(a) != get_num_clones(b)) {
        logIfVerbose(
            "number of clones unequal: ",
            "{get_num_clones(a)}, {get_num_clones(b)}"
        )
        return(FALSE)
    }

    if (!setequal(get_radii(a), get_radii(b))) {
        logIfVerbose("radii are not set-equal")
        logSetDiffIfVerbose(
            get_radii(a), get_radii(b),
            deparse(substitute(a)), deparse(substitute(b))
        )
        return(FALSE)
    }

    if (!setequal(get_clonotypes(a), get_clonotypes(b))) {
        logSetDiffIfVerbose(
            get_clonotypes(a), get_clonotypes(b),
            deparse(substitute(a)), deparse(substitute(b))
        )
        return(FALSE)
    }

    list_of_normed_ab <- list(a, b) %>%
        lapply(function(x)
            normalizeClusterListRadii(x, rad_decrease, clone_scale_factor)
        )
    
    list_of_normed_ab %>%
        append(list("threshold" = tolerance, "verbose" = verbose)) %>%
        applyListAsArgsTo(areGeometricallyEqualNormalizedClusterLists)

}

normalizeClusterListRadii <- function(
    clusterlist, rad_decrease, clone_scale_factor
) {
    set_radii(
        clusterlist,
        ((get_radii(clusterlist) + rad_decrease)^2) / clone_scale_factor
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
areGeometricallyEqualNormalizedClusterLists <- function(
    a,
    b,
    threshold = 1e-6,
    verbose = FALSE,
    left_name = "object",
    right_name = "expected"
) {

    logIfVerbose <- function(...) if (verbose) message(glue(...))

    list_of_clusterlist_dfs <- list(a, b) %>% lapply(function(clusterlist) {
        clusterlist %>%
            convert_to_dataframe("placeholder") %>%
            dplyr::mutate(r = r / min(r)) %>%
            dplyr::select(-label)
    })

    for (radius in sort(unique(get_radii(a)))) {
        
        logIfVerbose("checking for radius {radius}")
        
        clusterlist_dfs_filtered_by_curr_rad <- list_of_clusterlist_dfs %>%
            lapply(function(x) dplyr::filter(x, r == radius))

        are_dim_equal <- clusterlist_dfs_filtered_by_curr_rad %>%
            applyListAsArgsTo(function(a, b) identical(dim(a), dim(b)))

        if (!are_dim_equal) {
            logIfVerbose(
                "dimensions (clones, slots) unequal: ",
                "({nrow(a)}, {ncol(a)}), ({nrow(b), ncol(b)})"
            )
            return(FALSE)
        }

        clonotypes_list <- clusterlist_dfs_filtered_by_curr_rad %>%
                lapply(function(x) x[["clonotype"]])
        
        are_clonotypes_setequal <- clonotypes_list %>%
            applyListAsArgsTo(setequal)

        if (!are_clonotypes_setequal) {
            logIfVerbose(
                "clonotypes are not set-equal.\n\n",
                getSetDiffAsListInStr(
                    clonotypes_list[1], clonotypes_list[2],
                    left_name, right_name
                ),
            )
            return(FALSE)
        }

        xy_are_setequal <- clusterlist_dfs_filtered_by_curr_rad %>%
            lapply(function(a) dplyr::arrange(a, x, y)) %>%
            applyListAsArgsTo(function(df1, df2) {
                if (nrow(df1) + nrow(df2) == 0) return(TRUE)
                sum(df1 - df2) < (ncol(df1) * nrow(df1) * threshold)
            })

        if (!xy_are_setequal) {
            logIfVerbose("x, y coordinates are not set-equal")
            return(FALSE)
        }
    }

    return(TRUE)

}
