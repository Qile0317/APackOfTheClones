source("testdata/SCIP.R")
library(data.table)

#test_that("count_umap_clusters() works", {
#  expect_equal(count_umap_clusters(test_integrated_pbmc), 2)
#})

test_that("get_clone_sizes works", {
  trial <- get_clone_sizes(test_integrated_pbmc, scale_factor = 1)
  
  expect_equal(length(trial), 2)
  expect_equal(round(sum(trial[[1]])), 36)
  expect_equal(round(sum(trial[[2]])), 23)
  
  trial_cl_1_tabled <- as.numeric(unname(table(trial[[1]])))
  trial_cl_2_tabled <- as.numeric(unname(table(trial[[2]])))

  expect_equal(trial_cl_1_tabled, c(12,12,2,2))

  expect_equal(
    trial_cl_2_tabled, c(10,7,2))
})
