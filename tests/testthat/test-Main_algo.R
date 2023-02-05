source("testdata/cluster_lists.R")

# need to test prerequisite functions as well
# probably need another testdata file for nodes

test_that("est_rad() works", {
  expect_equal(est_rad(c1), c1[[5]])
  expect_equal(est_rad(c2), c2[[5]], tolerance = 3e-7)
  expect_equal(est_rad(c3), c3[[5]], tolerance = 5e-6)
})

# testing circle_layout
test_that("circle_layout() works", {
  expect_equal(circle_layout(c(1, 1, 1, 1, 1.3, 1, 1, 1, 1.3),
                             c1$centroid, progbar = FALSE),
               c1, tolerance = 3e-7)
  expect_equal(circle_layout(c(1.9, 1.5, 1.6, 1.5, 2.0, 1.5),
                             c2$centroid, progbar = FALSE),
               c2, tolerance = 3e-7)
  expect_false(identical(c3, circle_layout(c(0.8,0.6,0.6,0.5,0.5),
                                           progbar = FALSE)))
})

test_that("circle_layout() handles edge cases", {
  expect_null(circle_layout(list()))

  #input length = 1
  expect_equal(circle_layout(c(69)),
               list(
                 "x" = 0,
                 "y" = 0,
                 "rad" = 69,
                 "centroid" = c(0, 0),
                 "clRad" = 69))
  expect_equal(circle_layout(c(420), c(6, 9)),
               list(
                 "x" = 6,
                 "y" = 9,
                 "rad" = 420,
                 "centroid" = c(6, 9),
                 "clRad" = 420))

  # input length = 2
  expect_equal(circle_layout(c(69, 420)),
               list("x" = c(-420, 69),
                    "y" = c(0, 0),
                    "rad" = c(420, 69),
                    "centroid" = c(0, 0),
                    "clRad" = 244.5))
  expect_equal(circle_layout(c(69,420), c(420, 69)),
               list("x" = c(0, 489),
                    "y" = c(69, 69),
                    "rad" = c(420, 69),
                    "centroid" = c(420, 69),
                    "clRad" = 244.5))

  # input length = 3 (not really an edge case but rather a test of the while loop)
  expect_equal(circle_layout(c(1, 2, 3)),
               list("x" = c(-2.7333333, 2.2666667, 0.4666667),
                    "y" = c(-0.8, -0.8, 1.6),
                    "rad" = c(3, 2, 1),
                    "centroid" = c(0, 0),
                    "clRad" = 4.266667),
               tolerance = 3e-7)
})

# need to test the optional args

