test_that("find_centroid() works", {
  input <- load("tests/UMAP_coords.Rda")
  expected_output <- load("tests/expected_find_centroid_output.Rda")
  expect_equal(input,expected_output)
})

# unfinished
