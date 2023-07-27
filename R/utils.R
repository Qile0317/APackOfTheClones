# This dummy function definition is included with the package to ensure that
# 'tools::package_native_routine_registration_skeleton()' generates the required
# registration info for the 'run_testthat_tests' symbol
(function() {
    .Call("run_testthat_tests", FALSE, PACKAGE = "APackOfTheClones")
})

# progress bar functions

progress_bar <- function (x = 0, max = 100) {
    percent <- 100 * (x / max)
    cat(sprintf(
        '\r[%-50s] %d%%',
        paste(rep('=', percent * 0.5), collapse = ''),
        floor(percent)
    ))
}

start_progress_bar <- function(verbose = TRUE) {
    if (verbose) {
        progress_bar(0, 1)
    }
}

end_progress_bar <- function(verbose = TRUE) {
    if (verbose) {
        progress_bar(1, 1)
    }
}

print_completion_time <- function(start_time, digits = 2) {
    end_time <- Sys.time()
    message(paste(
        "\nCompleted successfully, time elapsed:",
        round(as.numeric(end_time - start_time), digits),
        "seconds\n"
    ))
}

# readability functions

isnt_empty <- function(inp) {
    !identical(inp, list())
}

isnt_na <- function(inp) {
    !any(is.na(inp))
}

isnt_empty_nor_na <- function(inp) {
    isnt_empty(inp) && isnt_na(inp)
}

is_int <- function(num) {
    return(num == as.integer(num))
}

should_estimate <- function(obj, auto_str = "auto") {
    identical(obj, auto_str)
}

get_xr <- function(plt) {
    ggplot2::ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range
}

get_yr <- function(plt) {
    ggplot2::ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range
}

# string related functions

attempt_correction <- function(s) {
    s <- tolower(s)
    if (identical(s, "t-sne")) {
        s <- "tsne"
    }
    s
}

closest_word <- function(s, strset = c("umap", "tsne", "pca")) {
    closest_w <- strset[1]
    closest_dist <- utils::adist(s, closest_w)
    for(i in 2:length(strset)) {
        curr_dist <- utils::adist(s, strset[i])
        if (curr_dist < closest_dist) {
            closest_w <- strset[i]
            closest_dist <- curr_dist
        }
    }
    closest_w
}

# utility
find_first_non_empty <- function(l) {
    for (item in l) {
        if (isnt_empty(item)) {
            return(item)
        }
    }
    return(NULL)
}

# R interface function for checking if metadata names to be added overlaps with
# existing (probably placed in wrong file :/)
metadata_name_warnstring <- function(seurat_obj, tcr_dataframe) {

    seurat_names <- names(seurat_obj@meta.data)
    tcr_names <- names(tcr_dataframe)

    if (any(is.na(tcr_names))) {
        return("tcr_dataframe has NAs in names, please fix")
    }

    if (has_common_strs(seurat_names, tcr_names)) {
        return("tcr_dataframe has repeated names with the seurat_object metadata")
    }

    return(NULL)
}
