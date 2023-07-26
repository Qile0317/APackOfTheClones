source("testdata/SCIP.R")

untable <- function(a) {
    as.numeric(unname(a))
}

test_that("add_raw_clone_sizes works", {
    test_apotc_obj <- initialize_apotc(2,1,1,'umap')
    expect_equal(
        add_raw_clone_sizes(
            initialize_apotc(2, 1, 1, 'umap'), test_integrated_pbmc
        ),
        test_apotc_obj_with_raw_clone_sizes
    )
})

test_that("get_transformed_clone_sizes works", {
    sizelist <- list(c(1, 1, 2), NULL, c(0.5, 3))
    expect_equal(
        get_transformed_clone_sizes(sizelist, 0.5, 3),
        list(c(0.5, 0.5, sqrt(2)*0.5), list(), c(sqrt(0.5)*0.5, sqrt(3)*0.5))
    )
})

# NEED TO

test_that("get_clone_sizes works", {
    trial <- get_clone_sizes(test_integrated_pbmc, scale_factor = 1)

    expect_equal(length(trial), 2)
    expect_equal(round(sum(trial[[1]])), 36)
    expect_equal(round(sum(trial[[2]])), 23)

    trial_cl_1_tabled <- untable(table(trial[[1]]))
    trial_cl_2_tabled <- untable(table(trial[[2]]))

    expect_equal(trial_cl_1_tabled, c(12, 12, 2, 2))
    expect_equal(trial_cl_2_tabled, c(10, 7, 2))
})

test_that("count_clone_sizes works", {
    trial <- count_clone_sizes(test_integrated_pbmc)
    expect_equal(length(trial), 2)
    expect_equal(untable(trial[[1]]), c(12, 12, 2, 2))
    expect_equal(untable(trial[[2]]), c(10, 7, 2))
})


