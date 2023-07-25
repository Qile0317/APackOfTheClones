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

process_rad_vec <- function(input_rad_vec, ORDER, scramble) {
    if (ORDER) {
        return(sort(input_rad_vec, decreasing = TRUE, method = "radix"))
    }
    if (scramble) {
        set.seed(42)
        return(sample(input_rad_vec))
    }
}

# vectorized circle_layout - outputs list of clusterlists
pack_into_clusterlists <- function(
    sizes, centroids, num_clusters, rad_scale = 1,
    ORDER = TRUE, scramble = FALSE, try_place = FALSE, verbose = TRUE
){
  output_list <- list()

  for(i in 1:num_clusters){
      input_rad_vec <- sizes[[i]]

      if (!isnt_empty(input_rad_vec) || is.null(input_rad_vec)) {
          output_list[[i]] <- list()
          next
      }

      if(verbose){
          message(paste("\npacking cluster", as.character(i-1)))
      }

      if (ORDER || scramble) { # assume only 1 is true
          input_rad_vec <- process_rad_vec(input_rad_vec, ORDER, scramble)
      }

      output_list[[i]] <- cpp_circle_layout(
          input_rad_vec,
          centroid = centroids[[i]],
          rad_scale_factor = rad_scale,
          try_place = try_place,
          verbose = verbose
      )
  }
  output_list
}
