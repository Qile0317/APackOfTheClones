data("mini_seurat_obj")

test_that("get_num_clusters works", {
  expect_equal(get_num_clusters(mini_seurat_obj), 2L)
})
