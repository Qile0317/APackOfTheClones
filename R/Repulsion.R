# Functions defined in src/repulsion.cpp
# - get_average_vector(vec_list)
#     avg vector of a list of vectors
# - get_component_repulsion_vector(inp,i,j,G):
#     repulsion vec of cluster i, j in inp
# - calculate_transformation_vectors(t,v,n):
#     get average vectors from a list of list of vecs
#     returns `list()` if all zeros
# - do_cluster_intersect(cn_c,cn_r,cm_c,cm_r,thr):
#     check if two clusterlists overlap in c++
# - do_proceed
#     check if two clusters need to be repulsed

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

# do_cl_intersect has been replaced by do_cluster_intersect, this is a wrapper 
do_cl_intersect <- function(Cn, Cm, thr = 1) {
  do_cluster_intersect(
    Cn[[4]], Cn[[5]], Cm[[4]], Cm[[5]], thr
  )
}

isnt_empty <- function(inp) !identical(inp, list())

do_proceed <- function(inp, i, j, thr) {
  if (i != j) {
    if (isnt_empty(inp[[i]])) {
      if (isnt_empty(inp[[j]])) {
        return(do_cl_intersect(inp[[i]], inp[[j]], thr))
      }
    }
  }
  return(FALSE) 
}

# O(N^2) operation to calculate all repulsion vectors for each cluster
calculate_repulsion_vectors <- function(
  overall_repulsion_vec, inp,
  num_clusters, G = 1, thr = 0
) {
  for (i in 1:num_clusters) { 
    for (j in 1:num_clusters) {
      if (do_proceed(inp,i,j,thr)) {
        overall_repulsion_vec[[i]][[j]] <- get_component_repulsion_vector(
          inp, i, j, G
        )
      }else {
        overall_repulsion_vec[[i]][[j]] <- c(0, 0)
      }
    }
  }
  overall_repulsion_vec
}

# iterative repulsion. inp is a list of clusterlists.
# missing members of clusterlists are NA right now
# near future: make it modify an apotc class and in Rcpp
repulse_cluster <- function(
  inp, thr = 1, G = 1, max_iter = 20, verbose = TRUE
) {
  if (G <= 0) {stop("repulsion strength must be a positive real number")}
  
  #init variables - could use a class
  num_clusters <- length(inp)
  transformation_vectors <- initialize_direction_vectors(num_clusters) # variable naming is confusing here; this is a list of the transformations for each cluster at the end of each iteration. 
  overall_repulsion_vec <- initialize_list_of_transformation_vectors(transformation_vectors, num_clusters) # this one is for storing all repulsion vectors for all pairwise comparisons that are yet to be averaged for each iteration

  for(curr_iteration in 1:max_iter){ 
    overall_repulsion_vec <- calculate_repulsion_vectors(
      overall_repulsion_vec, inp, num_clusters, G, thr
    )
    transformation_vectors <- calculate_transformation_vectors(
      transformation_vectors, overall_repulsion_vec, num_clusters
    )
    
    #transformation vectors is an empty list() if everything was c(0,0)
    if (identical(transformation_vectors, list())) {
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