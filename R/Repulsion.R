# Functions defined in src/repulsion.cpp
# pdV(c1,c2) | get the polar repulsion vector of 2 clusterlists
# comV(v) | convert from polar to component form
# get_average_vector(vec_list) | avg vector of a list of vectors

# compute component form of repulsion vector between two clusters in the
# clusterlist `inp`, assuming do_proceed(inp,i,j) == TRUE
# uses a modified version of coulombs law except its addition in the numerator
# to compute the magnitude of the repulsion vector / 2 with the same direction
get_component_repulsion_vector <- function(inp, i, j, G, dist_adjust = 0) {
  
  # find polar distance vector between centroids
  polar_dist_vec <- pdV(inp[[i]], inp[[j]])
  
  #Find polar repulsion vec - half the magnitude due to repeated i and j
  polar_repulsion_vec <- c(
    0.5 * G * (inp[[i]][[5]] + inp[[j]][[5]]) /
      (unname(polar_dist_vec["magnitude"]) + dist_adjust)^2,
    unname(polar_dist_vec["direction"])
  )
  
  # return the component form
  comV(polar_repulsion_vec)
}

# Alias to initialize direction vectors in a list
initialize_direction_vectors <- function(num_clusters) {
  direction_vectors <- vector("list", num_clusters)
  for (i in 1:num_clusters) {
    direction_vectors[[i]] <- c(0, 0)
  }
  direction_vectors
}

# Alias to initialize the overall repulsion vec
initialize_list_of_transformation_vectors <- function(blank_vectors, num_clusters) {
  output <- list()
  for (i in 1:num_clusters) {
    output[[i]] <- blank_vectors
  }
  output
}

#function to check if 2 cluster lists overlap, with a threshold.
do_cl_intersect <- function(Cn, Cm, thr = 1) {
  centroid_xdif <- (Cn$centroid[1] - Cm$centroid[1])
  centroid_ydif <- (Cn$centroid[2] - Cm$centroid[2])
  centroid_euc_dist <- sqrt((centroid_xdif^2) + (centroid_ydif^2))
  return(identical((centroid_euc_dist + thr) < (Cn$clRad + Cm$clRad), TRUE)) #idk why without identical() it returns logical(0) when false
}

do_proceed <- function(inp, i, j, thr) {
  if (i != j) {
    if (!any(is.na(inp[[i]]))) {
      if (!any(is.na(inp[[j]]))) {
        return(do_cl_intersect(inp[[i]], inp[[j]], thr))
      }
    }
  }
  return(FALSE) 
}

# O(N^2) operation to calculate all repulsion vectors for each cluster
calculate_repulsion_vectors <- function(overall_repulsion_vec, inp, num_clusters, G = 1, thr = 0, dist_adjust = 0) {
  for (i in 1:num_clusters) { 
    for (j in 1:num_clusters) {
      if (do_proceed(inp,i,j,thr)) {
        overall_repulsion_vec[[i]][[j]] <- get_component_repulsion_vector(
          inp, i, j, G, dist_adjust
        )
      }else {
        overall_repulsion_vec[[i]][[j]] <- c(0, 0)
      }
    }
  }
  overall_repulsion_vec
}

# function to just get the average vectors from a list of list of repulsion vectors
calculate_transformation_vectors <- function(transformation_vectors, overall_repulsion_vec, num_clusters) {
  contains_nonzero_vector <- FALSE
  for (i in 1:num_clusters) {
    transformation_vectors[[i]] <- get_average_vector(overall_repulsion_vec[[i]])
    if (!identical(transformation_vectors[[i]], c(0,0))) {
      contains_nonzero_vector <- TRUE
    }
  }
  if (contains_nonzero_vector) {
    return(transformation_vectors)
  }
  "no_change"
}

# iterative repulsion. inp is a list of clusterlists. works but still destroys the structure
repulse_cluster <- function(
  inp, thr = 1, G = 1, max_iter = 20, dist_adjust = 0, verbose = TRUE
  ) {
  if (G <= 0) {stop("repulsion strength must be a positive real number")}
  
  #init variables - could use a class
  num_clusters <- length(inp)
  transformation_vectors <- initialize_direction_vectors(num_clusters) # variable naming is confusing here; this is a list of the transformations for each cluster at the end of each iteration. 
  overall_repulsion_vec <- initialize_list_of_transformation_vectors(transformation_vectors, num_clusters) # this one is for storing all repulsion vectors for all pairwise comparisons that are yet to be averaged for each iteration

  for(curr_iteration in 1:max_iter){ 
    overall_repulsion_vec <- calculate_repulsion_vectors(
      overall_repulsion_vec, inp, num_clusters, G, thr, dist_adjust
    )
    transformation_vectors <- calculate_transformation_vectors(
      transformation_vectors, overall_repulsion_vec, num_clusters
    )
    if (identical(transformation_vectors, "no_change")) {
      if(verbose) {progress_bar(1,1)}
      return(inp)
    }
    # with the transformation vectors established, each cluster is moved
    for (i in 1:num_clusters) {
      if (!any(is.na(inp[[i]]))) {
        inp[[i]] <- trans_coord(inp[[i]], transformation_vectors[[i]])
      }
    }
    if (verbose) {
      progress_bar(curr_iteration, max_iter)
    }
  }
  inp
}

# could be rewritten from scratch in rust