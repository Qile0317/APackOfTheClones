source("testdata/cluster_lists.R")

# numerical operations
test_that("distV() can get direction vector from clusterlists", {
  expect_equal(distV(c1, c2), c(-3, -4))
})

test_that("polV() can do polar form conversion", {
  expect_equal(polV(c(1,1)), c("magnitude" = 1.41421356,
                               "direction" = 0.78539816))
})

test_that("pdV() can get polar direction vector", {
  expect_equal(pdV(c1, c2), c("magnitude" = 5,
                              "direction" = -2.2142974),
               tolerance = 3e-6)
  expect_equal(pdV(c1, c2), polV(distV(c1, c2)))
})

test_that("comV() can convert from polar to component form", {
  expect_equal(comV(c("magnitude" = 1.41421356,
                      "direction" = 0.78539816)),
               c(1, 1))
})

test_that("sumL() can sum component vectors in a list", {
  expect_equal(sumL(list(c(1, 1), c(1, 2), c(-1, -9))), c(1, -6))
})

# test cluster intersection function
test_that("do_cl_intersect() = TRUE for overlapping lists", {
  expect_true(do_cl_intersect(c1, c2))
  expect_true(do_cl_intersect(c1, c1_shifted_by_4_5))
})

test_that("do_cl_intersect() = FALSE for overlapping lists with high threshold", {
  expect_false(do_cl_intersect(c1, c2, 100))
})

test_that("do_cl_intersect() = FALSE for non-overlapping lists", {
  expect_false(do_cl_intersect(c1, c3))
  expect_false(do_cl_intersect(c2, c3))
  expect_false(do_cl_intersect(c1, c1_shifted_to_9_0))
})

test_that("do_cl_intersect() handles edge case", {
  expect_false(do_cl_intersect(NA, Cm))
})

# test the helper function for repulsion iteration
test_that("do_proceed() works", {
  expect_false(do_proceed(test_cluster_lists, 1, 1, 2))
  expect_true(do_proceed(test_cluster_lists, 1, 4, 2))
})

test_that("do_proceed() handles NAs", {
  new_list <- test_cluster_lists
  new_list[[3]] <- NA
  expect_false(do_proceed(new_list, 1, 3, 2))
  expect_true(do_proceed(new_list, 1, 4, 2))
  rm("new_list")
})
