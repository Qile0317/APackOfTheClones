# vectorized circle_layout - outputs list of clusterlists
pack_into_clusterlists <- function(
  sizes, centroids, num_clusters, rad_scale = 1,
  ORDER = TRUE, try_place = FALSE, verbose = TRUE
){
  output_list <- list()
  for(i in 1:num_clusters){
    input_rad_vec <- sizes[[i]]
    
    if (is.null(input_rad_vec) || identical(input_rad_vec, c(0))) {
      output_list[[i]] <- list()
    }else{
      
      if(verbose){
        message(paste("\npacking cluster", as.character(i-1)))
      }
      
      if (ORDER) {
        input_rad_vec <- sort(
          input_rad_vec, decreasing = TRUE, method = "radix"
        )
      }
      
      output_list[[i]] <- cpp_circle_layout( 
        input_rad_vec,
        centroid = centroids[[i]],
        rad_scale_factor = rad_scale,
        try_place = try_place,
        verbose = verbose
      )
    }
  }
  output_list
}
