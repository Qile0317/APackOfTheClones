sourcedata("v0", "cluster_lists")

# testing cpp_circle_layout
test_that("cpp_circle_layout() works", {
  skip_on_cran()

  expect_equal(cpp_circle_layout(c(1.3, 1.3, rep(1, 7)),
                             c1$centroid, verbose = FALSE),
               c1, tolerance = 1e-9)

  expect_equal(cpp_circle_layout(c(2.0, 1.9, 1.6, 1.5, 1.5, 1.5),
                             c2$centroid, verbose = FALSE),
               c2, tolerance = 1e-9)

  expect_false(identical(
    c3,
    cpp_circle_layout(
      c(0.8,0.6,0.6,0.5,0.5), centroid = c(0, 0), verbose = FALSE
    ),
  ))

  expect_equal(
    cpp_circle_layout(
      c(6,9,4,2,7), centroid = c(0,0),
      verbose = FALSE, try_place = FALSE
    ),

    list(
      x = c(
        -6.73333333333333, 8.26666666666667, -1.53333333333333,
        -7.50533138211336, -8.70501097444157
      ),
      y = c(
        2.84722086720835, 2.84722086720835, -5.6944417344167, -5.11544319783168,
        -14.0351275462686
      ),
      rad = c(6, 9, 4, 2, 7),
      centroid = c(0, 0),
      clRad = 17.2666666666667
    ),

    tolerance = 3e-6
  )
})

test_that("cpp_circle_layout() handles edge cases", {
  skip_on_cran()

  #input length = 1
  expect_equal(cpp_circle_layout(c(69), centroid = c(0,0), verbose = FALSE),
               list(
                 "x" = 0,
                 "y" = 0,
                 "rad" = 69,
                 "centroid" = c(0, 0),
                 "clRad" = 69)
               )
  expect_equal(cpp_circle_layout(c(420), c(6, 9), verbose = FALSE),
               list(
                 "x" = 6,
                 "y" = 9,
                 "rad" = 420,
                 "centroid" = c(6, 9),
                 "clRad" = 420)
               )

  # input length = 2
  expect_equal(cpp_circle_layout(c(69, 420), centroid = c(0,0), verbose = FALSE),
               list("x" = c(-69, 420),
                    "y" = c(0, 0),
                    "rad" = c(69, 420),
                    "centroid" = c(0, 0),
                    "clRad" = 244.5
                    )
               )
  expect_equal(cpp_circle_layout(c(420, 69), c(420, 69), verbose = FALSE),
               list("x" = c(0, 489),
                    "y" = c(69, 69),
                    "rad" = c(420, 69),
                    "centroid" = c(420, 69),
                    "clRad" = 244.5
                    )
               )

  # input length = 3
  expect_equal(cpp_circle_layout(c(3, 2, 1), centroid = c(0,0), verbose = FALSE),
               list(x = c(-2.73333333333333, 2.26666666666667, 0.466666666666667
               ), y = c(0.8, 0.8, -1.6), rad = c(3, 2, 1), centroid = c(0, 0
               ), clRad = 4.26666666666667),
               tolerance = 3e-7)
})

# Testing optional args (need more)
test_that("cpp_circle_layout(rad_decrease = 0.1) works", {
    skip_on_cran()

    new_c1 <- c1
    new_c1$rad <- c1$rad - 0.1
    new_c1$clRad <- c1$clRad - 0.1
    expect_equal(
        cpp_circle_layout(
            c(1.3, 1.3, rep(1, 7)),
            centroid = new_c1$centroid,
            rad_decrease = 0.1,
            verbose = FALSE
        ),
       new_c1,
       tolerance = 1e-8
    )
})

test_that("process_rad_vec works", {
    set.seed(42)
    expect_equal(process_rad_vec(c(4,2,5,1,3), TRUE, FALSE), c(5,4,3,2,1))
    expect_equal(process_rad_vec(c(4,2,5,1,3), FALSE, TRUE), c(4,3,1,5,2))
    expect_equal(process_rad_vec(c(4,2,5,1,3), FALSE, FALSE), c(4,2,5,1,3))
})

test_that("pack_into_clusterlists works", {
    skip_on_cran()

    expect_equal(
        pack_into_clusterlists(
            sizes = test_radii,
            centroids = test_centroids,
            num_clusters = 5,
            verbose = FALSE
        ),
        test_cluster_lists,
        tolerance = 1e-8
    )

    expect_equal(
        pack_into_clusterlists(
            test_radii[1:4],
            test_centroids[1:4],
            4,
            verbose = FALSE
        ),
        test_cluster_lists[1:4],
        tolerance = 1e-8
    )

    expect_equal(
        pack_into_clusterlists(
            test_radii[1:4],
            test_centroids[1:4],
            4,
            verbose = FALSE
        ),
        pack_into_clusterlists(
            test_radii,
            test_centroids,
            5,
            verbose = FALSE
        )[1:4],
        tolerance = 1e-8
    )
})

test_that("pack_into_clusterlists handles NULLS", {
    skip_on_cran()
    
    expect_equal(
        pack_into_clusterlists(
            list(NULL, c1[[3]]),
            list(c(0, 0), c(0, 0)),
            2,
            verbose = FALSE
        ),
        list(list(), c1),
        tolerance = 1e-8
    )
})
