data("mini_seurat_obj")

test_that("get_num_clusters works", {
  expect_equal(get_num_clusters(mini_seurat_obj), 2L)
})

test_that("isnt_empty works", {
  expect_true(isnt_empty(list(c(1,2,3))))
  expect_false(isnt_empty(list()))
})

test_that("isnt_na works", {
  expect_true(isnt_na(list(c(1,2,3))))
  expect_false(isnt_na(list(NA,NA)))
})

test_that("isnt_empty_nor_na works", {
  expect_true(isnt_empty_nor_na(list(c(1,2,3))))
  expect_false(isnt_empty_nor_na(list(NA,NA)))
  expect_false(isnt_empty_nor_na(list()))
})