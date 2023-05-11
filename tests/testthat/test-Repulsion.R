source("testdata/cluster_lists.R")

test_that("get_average_vector() works", {
  expect_equal(get_average_vector(list(c(1, 1), c(1, 2), c(-1, -9))), c(1/3, -2))
  expect_equal(get_average_vector(list(c(0,0),c(0,0),c(0,0))), c(0,0))
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

# need to test do_proceed

test_that("get_component_repulsion_vector works", {
  expect_equal(
    get_component_repulsion_vector(
      list(c1,c2), 1, 2, 1
    ),
    c(-0.1091579, -0.1455439),
    tolerance = 1e-5
  )
})

test_that("initialize_direction_vectors() works", {
  expect_equal(initialize_direction_vectors(2), list(c(0,0),c(0,0)))
})

test_that("initialize_list_of_transformation_vectors() works", {
  blank_vecs <- initialize_direction_vectors(2)
  expect_equal(
    initialize_list_of_transformation_vectors(blank_vecs,2),
    list(blank_vecs,blank_vecs))
})

# need to test more edgecases!
test_that("calculate_repulsion_vectors() works", {
  overall_repulsion_vec <- initialize_list_of_transformation_vectors(
    initialize_direction_vectors(3),3
  )
  inp <- list(c1, c1_shifted_by_4_5, c2)
  expect_equal(
    calculate_repulsion_vectors(
      overall_repulsion_vec, inp, 3
    ),
    list(
      list(c(0,0), c(-0.06768692, -0.08460865), c(-0.1091579, -0.1455439)),
      list(c(0.06768692, 0.08460865), c(0,0), c(1.608048, 1.608048)),
      list(c(0.1091579, 0.1455439), c(-1.608048, -1.608048), c(0,0))
    ),
    tolerance = 1e-6
  )
})

test_that("calculate_transformation_vectors() works", {
  expect_equal(
    calculate_transformation_vectors(
      initialize_direction_vectors(3),
      list(
        list(c(0,0), c(-0.1503472, -0.1879340), c(-0.2481042, -0.3308055)),
        list(c(0.1503472, 0.1879340), c(0,0), c(3.654919, 3.654919)),
        list(c(0.2481042, 0.3308055), c(-3.654919, -3.654919), c(0,0))
      ),
      3
    ),
    list(
      c(-0.1992257, -0.2593697),
      c(1.902633, 1.921426),
      c(-1.703407, -1.662057)
    ),
    tolerance = 1e-6
  )
})

test_that("repulse_cluster() works", {
  expect_equal(
    repulse_cluster(
      list(c1, c1_shifted_by_4_5),
      thr = 1, G = 1, max_iter = 10,
      verbose = F
    ),
    list(list(x = c(-1.792939452373, 0.807060547627004, -0.492939452372996, 
                    -2.47021185105545, 1.48433294630945, -0.492939452372996, 0.235321256461886, 
                    -1.74195114222057, 2.94949010969326), y = c(0.0162812165674311, 
                                                                0.0162812165674311, -1.88108537953359, -2.18174109860204, -2.18174109860204, 
                                                                1.91364781266846, -3.74378136652861, -4.04443708559705, -0.820376163478056
                    ), rad = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1), centroid = c(-0.492939452372996, 
                                                                            -0.616174315466245), clRad = 4.44242956206626), list(x = c(3.192939452373, 
                                                                                                                                       5.792939452373, 4.492939452373, 2.51566705369055, 6.47021185105545, 
                                                                                                                                       4.492939452373, 5.22120016120788, 3.24392776252543, 7.93536901443926
                                                                            ), y = c(6.24862984749993, 6.24862984749993, 4.3512632513989, 
                                                                                     4.05060753233045, 4.05060753233045, 8.14599644360095, 2.48856726440389, 
                                                                                     2.18791154533545, 5.41197246745444), rad = c(1.3, 1.3, 1, 1, 
                                                                                                                                  1, 1, 1, 1, 1), centroid = c(4.492939452373, 5.61617431546625
                                                                                                                                  ), clRad = 4.44242956206626)),
    tolerance = 1e-6 # dput() was used
  )
})