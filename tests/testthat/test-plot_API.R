source("testdata/cluster_lists.R")
suppressPackageStartupMessages(library(dplyr))
library(vdiffr)
library(ggforce)

test_that("df_full_join() works", {
  
  test_df <- structure(list(label = c("cluster 0", "cluster 0", "cluster 0", 
                                      "cluster 0", "cluster 0", "cluster 0", "cluster 0", "cluster 0", 
                                      "cluster 0", "cluster 1", "cluster 1", "cluster 1", "cluster 1", 
                                      "cluster 1", "cluster 1", "cluster 1", "cluster 1", "cluster 1", 
                                      "cluster 2", "cluster 2", "cluster 2", "cluster 2", "cluster 2", 
                                      "cluster 2", "cluster 2", "cluster 2", "cluster 2", "cluster 3", 
                                      "cluster 3", "cluster 3", "cluster 3", "cluster 3", "cluster 3", 
                                      "cluster 4", "cluster 4", "cluster 4", "cluster 4", "cluster 4"
  ), x = c(-1.3, 1.3, -1.48029736616688e-16, -1.97727239868245, 
           1.97727239868245, 2.22044604925031e-16, 0.728260708834882, -1.24901168984757, 
           3.44242956206626, 2.7, 5.3, 4, 2.02272760131755, 5.97727239868245, 
           4, 4.72826070883488, 2.75098831015243, 7.44242956206626, 7.7, 
           10.3, 9, 7.02272760131755, 10.9772723986824, 9, 9.72826070883488, 
           7.75098831015243, 12.4424295620663, 1.01965811965812, 4.91965811965812, 
           3.06068376068376, -0.0162220916735012, 6.15406282750692, 3.05811965811966, 
           19.2380952380952, 20.6380952380952, 20.1238095238095, 19.042322998886, 
           20.1095238095238), y = c(0.632455532033676, 0.632455532033676, 
                                    -1.26491106406735, -1.56556678313579, -1.56556678313579, 2.5298221281347, 
                                    -3.12760705106236, -3.4282627701308, -0.204201848011811, 5.63245553203368, 
                                    5.63245553203368, 3.73508893593265, 3.43443321686421, 3.43443321686421, 
                                    7.5298221281347, 1.87239294893764, 1.5717372298692, 4.79579815198819, 
                                    0.632455532033676, 0.632455532033676, -1.26491106406735, -1.56556678313579, 
                                    -1.56556678313579, 2.5298221281347, -3.12760705106236, -3.4282627701308, 
                                    -0.204201848011811, 4.98850135349066, 4.98850135349066, 2.02299729301867, 
                                    1.64530644625424, 1.82049803637739, 7.83361550870711, 20.3614031611621, 
                                    20.3614031611621, 19.2771936776758, 19.076228764298, 21.3260856192269
           ), r = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1, 1.3, 1.3, 1, 1, 1, 1, 
                    1, 1, 1, 1.3, 1.3, 1, 1, 1, 1, 1, 1, 1, 2, 1.9, 1.6, 1.5, 1.5, 
                    1.5, 0.8, 0.6, 0.6, 0.5, 0.5)), class = "data.frame", row.names = c(NA, 
                                                                                        -38L))

  expect_equal(
    df_full_join(
      list(c1,c1_shifted_by_4_5, c1_shifted_to_9_0, c2,c3)),
    test_df)

  expect_equal(df_full_join(list(c1)),
               test_df[1:9, ])
  
  expect_equal(df_full_join(list(c1, NULL)), test_df[1:9, ])
})

test_that("plot_clusters() works", {
  test_plot <- plot_clusters(
    insert_colors(df_full_join(list(c1,c2,c3)),3))
  expect_doppelganger("plot_c1_c2_c3", test_plot)
})

test_that("plot_API() works", {
  test_plot <- plot_API(
    test_radii[-1], test_centroids[-1], 5,
    progbar = F, repulse = F
  )
  expect_doppelganger("plot_API_all", test_plot)
})

test_that("plot_API(rad_decrease = 0.8) works", {
  test_plot <- plot_API(
    test_radii[-1], test_centroids[-1], 5,
    progbar = F, repulse = F,
    rad_decrease = 0.8
  )
  expect_doppelganger(
    "plot_API_all_scaled_0.8",
    test_plot
  )
})

test_that("plot_API(try_place = TRUE) works", {
  test_plot <- plot_API(
    test_radii[-1], test_centroids[-1], 5,
    progbar = F, repulse = F,
    try_place = T
  )
  expect_doppelganger(
    "plot_API_try_place_TRUE",
    test_plot
  )
})

# need to test more edge cases
