#progress_bar_with_text <- function(verbose, x, max, info, num = -1) {
#  percent <- 100 * (x / max)
#  if (num != -1) {
#    info <- paste(info, num)
#  }
#  cat(sprintf(
#    '\r%s [%-50s] %d%%',
#    info,
#    paste(rep('=', percent * 0.5), collapse = ''),
#    round(percent)
#  ))
#}

# vectorized circle_layout - outputs list of clusterlists
pack_into_clusterlists <- function(
  sizes, centroids, num_clusters, rad_scale = 1,
  ORDER = TRUE, scramble = FALSE, try_place = FALSE, verbose = TRUE
){
  #if (verbose) {
  #  progress_bar_with_text(0,1, "Initializing...")
  #}
  output_list <- list()
  for(i in 1:num_clusters){
    #if (verbose) {
    #  progress_bar_with_text(i, num_clusters, "Packing Cluster", i)
    #}
    
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
      if (scramble) {
        set.seed(42)
        input_rad_vec <- sample(input_rad_vec)
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
  #if (verbose) {
  #  progress_bar_with_text(1,1, "Packing complete")
  #}
  output_list
}
