# test_that() is the most basic testing func
# describe() should verify that the right things were implemented

# test data, three cluster lists. c1 and c2 overlap while c3 is far away
c1 <- list("x" = c(2.7, 5.3, 4, 4, 4, 5.977272, 2.022728, 3.271739, 5.249012),
           "y" = c(4.367544, 4.367544, 6.264911, 2.470178, 6.264911, 6.565567,
                   6.565567, 8.127607, 8.428263),
           "rad" = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1),
           "centroid" = c(4, 5),
           "clRad" = 2.977272)

c2 <- list("x" = c(1.0196581, 4.9196581, 3.0606838, 3.0581197, 3.0581197,
                   0.1059401),
           "y" = c(3.0114986, 3.0114986, 5.9770027, 0.1663845, 5.8566128,
                   6.3901255),
           "rad" = c(2, 1.9, 1.6, 1.5, 1.5, 1.5),
           "centroid" = c(3, 4),
           "clRad" = 3.819658)

c3 <- list("x" = c(19.2381, 20.6381, 20.12381, 20.10952, 20.10952),
           "y" = c(19.6386, 19.6386, 20.72281, 18.67391, 20.60328),
           "rad" = c(0.8, 0.6, 0.6, 0.5, 0.5),
           "centroid" = c(20, 20),
           "clRad" = 1.238095)

# plot_clusters(df_full_join(list(c1,c2,c3)))

# Vector operations - tbh they are simple enough that they don't really need testing
test_that("distV() can get direction vector from clusterlists", {
  expect_equal(distV(c1, c2), c(1, 1))
})

test_that("polV() can do polar form conversion", {
  expect_equal(polV(c(1,1)), c("magnitude" = 1.41421356,
                               "direction" = 0.78539816))
})

test_that("pdV() can get polar direction vector", {
  expect_equal(pdV(c1, c2), c("magnitude" = 1.41421356,
                              "direction" = 0.78539816))
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
})

test_that("do_cl_intersect() = FALSE for overlapping lists with high threshold", {
  expect_false(do_cl_intersect(c1, c2, 100))
})

test_that("do_cl_intersect() = FALSE for non-overlapping lists", {
  expect_false(do_cl_intersect(c1, c3))
  expect_false(do_cl_intersect(c2, c3))
})

# test the actual repulsion

