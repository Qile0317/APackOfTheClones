source("testdata/cluster_lists.R")

# need to test prerequisite functions as well
# probably need another testdata file for nodes

# remember it assumes centroid at origin
test_that("est_rad() works", {
  expect_equal(est_rad(c1), c1[[5]])
  expect_equal(est_rad(c2), c2[[5]], tolerance = 3e-7)
  expect_equal(est_rad(c3), c3[[5]], tolerance = 5e-6)
})

test_that("circle_layout() works", {
  expect_null(circle_layout(list()))
  # inplen = 1
  # inplen = 2
  # inplen = 3
  expect_equal(circle_layout(c(1, 1, 1, 1, 1.3, 1, 1, 1, 1.3),
                             c1$centroid), c1, tolerance = 3e-7)
  expect_equal(circle_layout(c(1.9, 1.5, 1.6, 1.5, 2.0, 1.5),
                             c2$centroid),
               c2, tolerance = 3e-7)
  expect_false(identical(c3, circle_layout(c(0.8,0.6,0.6,0.5,0.5))))
})

