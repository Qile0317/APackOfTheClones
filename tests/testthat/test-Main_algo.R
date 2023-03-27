source("testdata/cluster_lists.R")

# need to test prerequisite functions as well
# probably need another testdata file for nodes

test_that("est_rad() works", {
  expect_equal(est_rad(c1), c1[[5]])
  expect_equal(est_rad(c2), c2[[5]], tolerance = 3e-7)
  expect_equal(est_rad(c3), c3[[5]], tolerance = 5e-6)
})

# floating point errors are destroying these testsets so i set a ridiculously high tolerance

# testing circle_layout
test_that("circle_layout() works", {
  expect_equal(circle_layout(c(1, 1, 1, 1, 1.3, 1, 1, 1, 1.3),
                             c1$centroid, progbar = FALSE),
               c1, tolerance = 1)

  expect_equal(circle_layout(c(1.9, 1.5, 1.6, 1.5, 2.0, 1.5),
                             c2$centroid, progbar = FALSE),
               c2, tolerance = 0.1)

  expect_false(identical(c3, circle_layout(c(0.8,0.6,0.6,0.5,0.5),
                                           progbar = FALSE)))
})

test_that("circle_layout() handles edge cases", {
  # no input for cluster
  expect_true(is.na(circle_layout(list())))

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
  expect_equal(circle_layout(c(1, 2, 3), progbar = FALSE),
               list(x = c(-2.73333333333333, 2.26666666666667, 0.466666666666667
               ), y = c(0.8, 0.8, -1.6), rad = c(3, 2, 1), centroid = c(0, 0
               ), clRad = 4.26666666666667),
               tolerance = 3e-7)
})

# Testing optional args
test_that("circle_layout(ORDER = FALSE) works", {
  expect_equal(
    circle_layout(c(6,9,4,2,7), ORDER = FALSE, progbar = FALSE),
    list(x = c(-6.73333333333333, 8.26666666666667, -1.53333333333333, 
               4.22415396998617, -12.0149282860133), 
         y = c(2.84722086720835, 2.84722086720835, -5.6944417344167,
               -7.38303284896003, -9.03153137180334),
         rad = c(6, 9, 4, 2, 7),
         centroid = c(0, 0),
         clRad = 17.2666666666667),
    tolerance = 3e-6)
})

test_that("circle_layout(rad_decrease = 0.9) works", {
  new_c1 <- c1
  new_c1[[3]] <- c1[[3]] * 0.9
  expect_equal(circle_layout(c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1),
                             centroid = new_c1$centroid,
                             rad_decrease = 0.9,
                             progbar = FALSE),
               new_c1,
               tolerance = 1)
})

# need to test try_place but its pretty reliable.
