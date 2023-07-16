# vectorized circle_layout - outputs list of clusterlists
pack_into_clusterlists <- function(
  sizes, centroids, num_clusters, rad_scale = 1,
  ORDER = TRUE, try_place = FALSE, verbose = TRUE
){
  output_list <- list()
  for(i in 1:num_clusters){
    if (is.null(sizes[[i]]) || identical(sizes[[i]], c(0))) {
      output_list[[i]] <- list()
    }else{
      if(verbose){
        message(paste("\npacking cluster", as.character(i-1)))
      }
      output_list[[i]] <- cpp_circle_layout( 
        sizes[[i]],
        centroid = centroids[[i]],
        rad_scale_factor = rad_scale,
        ORDER = ORDER,
        try_place = try_place,
        verbose = verbose
      )
    }
  }
  output_list
}
