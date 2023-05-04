# it creates the OPPOSITE vector so there is no need for G to be negative
distV <- function(c1, c2){
  return(c(c1$centroid[1] - c2$centroid[1],
           c1$centroid[2] - c2$centroid[2]))
} # works

#polar form conversion from component form. its with respect to x axis. - works
polV <- function(vec) {
  return(c("magnitude" = sqrt(sum(vec^2)),
           "direction" = atan2(vec[2], vec[1])))
  }

#polar distance vector - works
pdV <- function(c1, c2) {
  polV(distV(c1, c2))
}

#converts polar form to component form of vector - works
comV <- function(Pvec) {
  unname(c(Pvec[1] * cos(Pvec[2]), Pvec[1] * sin(Pvec[2])))
}

# get the average vector (basically decrease magnitude geometrically) (inefficient)
get_average_vector <- function(vec_list) {
  sum_vector <- c(0, 0)
  num_non_zero_vectors <- 0
  for (vec in vec_list) {
    if (!identical(vec, c(0,0))) {
      num_non_zero_vectors <- num_non_zero_vectors + 1
      sum_vector <- sum_vector + vec
    }
  }
  if (num_non_zero_vectors > 0) {
    return(sum_vector / num_non_zero_vectors)
  }
  c(0, 0)
}

get_component_repulsion_vector <- function(inp, i, j, G, dist_adjust = 0) {
  # find polar distance vector between centroids
  polar_dist_vec <- pdV(inp[[i]], inp[[j]])
  
  #Find polar repulsion vec - half the magnitude due to repeated i and j
  polar_repulsion_vec <- c(
    0.5 * (G * (inp[[i]][[5]] + inp[[j]][[5]])) /
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
    direction_vectors[[i]]<-c(0, 0)
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
  # super dumb, temp fix
  if (any(is.na(Cn)) || any(is.na(Cm))) {
    return(FALSE)
  }
  if (any(is.null(Cn)) || any(is.null(Cm))) {
    return(FALSE)
  }

  #calculate euclidean distance of centroids
  centroid_xdif <- (Cn$centroid[1] - Cm$centroid[1])
  centroid_ydif <- (Cn$centroid[2] - Cm$centroid[2])
  centroid_euc_dist <- sqrt((centroid_xdif^2) + (centroid_ydif^2))
  
  # return
  return((centroid_euc_dist + thr) < (Cn$clRad + Cm$clRad))
}

do_proceed <- function(inp, i, j, thr) {
  if (i != j) {
    if (length(inp) != 0) {
      if (length(inp[[i]]) != 0) {
        if (length(inp[[j]]) != 0) {
          if (length(thr) != 0) {
            if (do_cl_intersect(inp[[i]], inp[[j]], thr)) {
              return(TRUE)
            }
          }
        }
      }
    }
  }
  FALSE
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
  if (G <= 0) {
    stop("repulsion strength must be a positive real number")
  }
  
  #init variables
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
      inp[[i]] <- trans_coord(inp[[i]],transformation_vectors[[i]])
    }
    if (verbose) {
      progress_bar(curr_iteration, max_iter)
    }
  }
  inp
}

# could be rewritten from scratch in rust - actually pretty fast lol
# there could be a stackoverflow for big inputs? a bajilion functions are repeated called :P