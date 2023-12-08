# vectorized circle_layout - outputs list of clusterlists
pack_into_clusterlists <- function(
    sizes, centroids, num_clusters, rad_decrease = 0,
    ORDER = TRUE, scramble = FALSE, try_place = FALSE, verbose = TRUE
){
  output_list <- vector("list", num_clusters)

  for(i in 1:num_clusters){
      input_rad_vec <- sizes[[i]]

      if (!isnt_empty(input_rad_vec) || is.null(input_rad_vec)) {
          output_list[[i]] <- list()
          next
      }

      if(verbose){message(paste("\npacking cluster", as.character(i-1)))}

      output_list[[i]] <- cpp_circle_layout(
          input_rad_vec = process_rad_vec(input_rad_vec, ORDER, scramble),
          centroid = centroids[[i]],
          rad_decrease = rad_decrease,
          try_place = try_place,
          verbose = verbose
      )
  }
  output_list
}

process_rad_vec <- function(input_rad_vec, ORDER, scramble) {
  if (ORDER) {
    return(sort(input_rad_vec, decreasing = TRUE))
  }
  if (scramble) {
    set.seed(42)
    return(sample(input_rad_vec))
  }
  input_rad_vec
}
