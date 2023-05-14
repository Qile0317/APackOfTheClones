source("testdata/cluster_lists.R")

test_that("get_average_vector() works", {
  expect_equal(get_average_vector(list(c(1, 1), c(1, 2), c(-1, -9))), c(1/3,-2))
  expect_equal(get_average_vector(list(c(0, 0), c(0, 0), c(0, 0))), c(0, 0))
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

test_that("do_proceed() works", {
  expect_false(do_proceed(test_cluster_lists, 1, 5, 0))
  expect_true(do_proceed(test_cluster_lists, 1, 2, 0))
  expect_false(do_proceed(test_cluster_lists, 1, 2, 20))
})

test_that("get_component_repulsion_vector works", {
  expect_equal(
    get_component_repulsion_vector(
      list(c1, c2), 1, 2, 1
    ),
    c(-0.1008981, -0.1345308),
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
    list(list(c(0, 0), c(-0.069912727392497, -0.0873909092406212), 
              c(-0.100898062586345, -0.134530750115127)), list(c(0.0699127273924969, 
                                                                 0.0873909092406212), c(0, 0), c(1.48636883882061, 1.48636883882061
                                                                 )), list(c(0.100898062586345, 0.134530750115127), c(-1.48636883882061, 
                                                                                                                     -1.48636883882061), c(0, 0))),
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
  
  expect_identical(
    calculate_transformation_vectors(
      initialize_direction_vectors(2),
      list(
        list(c(0,0),c(0,0)),
        list(c(0,0),c(0,0))
      ),
      2
    ),
    list()
  )
})

test_that("repulse_cluster() works", {
  expect_equal(
    repulse_cluster(
      list(c1, c1_shifted_by_4_5),
      thr = 1, G = 1, max_iter = 10,
      verbose = F
    ),
    list(list(x = c(-1.85085353933408, 0.749146460665919, -0.550853539334081, 
                    -2.52812593801653, -1.27911424816896, 0.698158150513489, 1.99266321363309, 
                    -0.550853539334081, 3.03766022320323), y = c(-0.0561113921339256, 
                                                                 -0.0561113921339256, -1.95347798823495, -2.25413370730339, -3.81617397522996, 
                                                                 -3.51551825616153, -1.99096696428653, 1.8412552039671, -0.285686562126976
                    ), rad = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1), centroid = c(-0.550853539334081, 
                                                                            -0.688566924167602), clRad = 4.58851376253731), list(x = c(3.25085353933408, 
                                                                                                                                       5.85085353933408, 4.55085353933408, 2.57358114065163, 3.8225928304992, 
                                                                                                                                       5.79986522918165, 7.09437029230125, 4.55085353933408, 8.13936730187139
                                                                            ), y = c(6.32102245620128, 6.32102245620128, 4.42365586010025, 
                                                                                     4.12300014103181, 2.56095987310524, 2.86161559217367, 4.38616688404867, 
                                                                                     8.2183890523023, 6.09144728620823), rad = c(1.3, 1.3, 1, 1, 1, 
                                                                                                                                 1, 1, 1, 1), centroid = c(4.55085353933408, 5.6885669241676), 
                                                                            clRad = 4.58851376253731)),
    tolerance = 1e-6 # dput() was used
  )
})