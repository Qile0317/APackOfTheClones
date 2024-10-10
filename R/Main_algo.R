# vectorized circle_layout - outputs list of clusterlists
pack_into_clusterlists <- function(
    sizes, centroids, num_clusters, rad_decrease = 0,
    ORDER = TRUE, scramble = FALSE, try_place = FALSE, verbose = TRUE
){

    output_list <- init_list(num_clusters, init_val = list())

    if (verbose) {
        packed_clone_count <- 0
        total_clone_count <- sum(sapply(sizes, length))
        start_progress_bar()
    }
  
    for (i in 1:num_clusters){
        input_rad_vec <- sizes[[i]]

        if (is_empty(input_rad_vec)) next

        processed_rad_vec <- process_rad_vec(input_rad_vec, ORDER, scramble)

        output_list[[i]] <- cpp_circle_layout(
            input_rad_vec = processed_rad_vec,
            centroid = centroids[[i]],
            rad_decrease = rad_decrease,
            try_place = try_place,
            verbose = FALSE
        )

        output_list[[i]][["clonotype"]] <- names(processed_rad_vec)

        if (verbose) {
          packed_clone_count %+=% length(input_rad_vec) # nolint
          progress_bar(packed_clone_count, total_clone_count)
        }
    }

    if (verbose) message("")
    output_list
}

process_rad_vec <- function(input_rad_vec, ORDER, scramble) {
    if (ORDER) {
        return(sort(input_rad_vec, decreasing = TRUE))
    }
    if (scramble) {
        return(sample(input_rad_vec))
    }
    input_rad_vec
}
