data("mini_seurat_obj")

test_that("progress_bar works", {
  expect_identical(
    capture.output(progress_bar(1, 1)),
    "\r[==================================================] 100%"
  )
  expect_identical(
    capture.output(progress_bar(69, 420)),
    "\r[========                                          ] 16%"
  )
})

test_that("isnt_empty works", {
  expect_true(isnt_empty(list(c(1,2,3))))
  expect_false(isnt_empty(list()))
})

test_that("isnt_na works", {
  expect_true(isnt_na(list(c(1, 2, 3))))
  expect_false(isnt_na(list(NA, NA)))
})

test_that("isnt_empty_nor_na works", {
  expect_true(isnt_empty_nor_na(list(c(1, 2, 3))))
  expect_false(isnt_empty_nor_na(list(NA, NA)))
  expect_false(isnt_empty_nor_na(list()))
})

test_that("attempt_correction works", {
  expect_identical("umap", attempt_correction("Umap"))
  expect_identical("tsne", attempt_correction("t-SNE"))
  expect_identical("pca", attempt_correction("PCA"))
})