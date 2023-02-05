suppressPackageStartupMessages(library(hash))

test_that("dict() syntax works", {
  expected_output <- hash()
  expected_output[["test1"]] <- "success"
  expected_output[["test2"]] <- "success"
  expect_equal(dict(c("test1", "test2"),c("success","success")), expected_output)
})

# cant test alot of the other functions without a seurat object...

test_that("dict_list_vals() works", {
  test_dict <- dict(c("test1", "test2"), c(2, 4))
  testing_dict_list <- list(test_dict)[rep(1,5)]
  expected_output <- list(c(0.002, 0.004))[rep(1,5)]

  expect_equal(dict_list_vals(testing_dict_list),expected_output)
})
