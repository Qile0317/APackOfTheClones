source("testdata/SCIP.R")

test_that("count_umap_clusters() works", {
  expect_equal(count_umap_clusters(test_integrated_pbmc), 2)
})

test_that("get_clone_sizes works", {
  expect_identical(
    get_clone_sizes(test_integrated_pbmc, scale_factor = 1),
    list(
      c(1, 1, 2, 3, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 4, 4, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 3),
      c(1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 2, 1, 1, 1, 2, 3, 2, 2, 2)
    )
  )

  expect_equal(
    get_clone_sizes(test_integrated_pbmc, scale_factor = 0.1),
    list(
      c(0.1, 0.1, 0.2, 0.3, 0.1, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.2, 0.1, 0.2,
        0.4, 0.4, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.2, 0.2, 0.3),
      c(0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2, 0.3, 0.2, 0.2, 0.2)
    )
  )
})

new_get_clone_sizes(test_integrated_pbmc)
