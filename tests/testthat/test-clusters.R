sourcedata("v0", "cluster_lists")
sourcedata("v0", "SCIP")

# TODO find_centroids

test_that("get_cluster_centroids() works", {
  expect_equal(
    get_cluster_centroids(test_pbmc),
    list(
      c(5.72926706448197, 8.48345941543579),
      c(-9.54877844080329, -14.1390990257263)
    ),
    tolerance = 1e-9
  )
  
  expect_equal(
    get_cluster_centroids(test_pbmc, "tsne"),
    list(
      c(-7.762187, 13.409124),
      c(12.93698, -22.34854)
    ),
    tolerance = 1e-6
  )
})

test_that("trans_coord() works withOUT new_coord", {
  c1_new <- c1
  c1_new$centroid <- c(4, 5)
  expect_equal(
    trans_coord(c1_new), c1_shifted_by_4_5, tolerance = 1e-9
  )
})

test_that("trans_coord() works WITH new_coord", {
  expect_equal(trans_coord(c1, c(9, 0)), c1_shifted_to_9_0, tolerance = 1e-9)
})

# need to test if it doesn't work when the input is null/na!

test_that("move_cluster() works", {
  expect_equal(move_cluster(c1, c(4, 5)), c1_shifted_by_4_5, tolerance = 1e-9)
  expect_equal(move_cluster(c1, c(9, 0)), c1_shifted_to_9_0, tolerance = 1e-9)
})
