source("testdata/cluster_lists.R")

test_that("df_full_join() works", {
  test_df <- readRDS("testdata/df_full_join_test_df.rds")
    
  expect_equal(
    df_full_join(
      list(c1,c1_shifted_by_4_5, c1_shifted_to_9_0, c2,c3)),
    test_df)

  expect_equal(df_full_join(list(c1)),
               test_df[1:9, ])
  
  expect_equal(df_full_join(list(c1, NA)), test_df[1:9, ])
  
  test_df_rows_1_till_9_with_label_1 <- test_df[1:9,]
  test_df_rows_1_till_9_with_label_1[[1]] <- rep("cluster 1", 9)
  
  expect_equal(
    df_full_join(list(NA, c1)),
    test_df_rows_1_till_9_with_label_1,
    tolerance = 1e-6
  )
  
  expect_identical(
    test_df_rows_1_till_9_with_label_1,
    df_full_join(list(NA, c1, NA))
  )
  
  expect_identical(
    df_full_join(list(NA, c1, NA, c2, NA)),
    df_full_join(list(list(), c1, list(), c2, list()))
  )
})

test_that("plot_clusters() works", {
  test_plot <- plot_clusters(
    insert_colors(df_full_join(list(c1,c2,c3)),3))
  expect_doppelganger("plot_c1_c2_c3", test_plot)
  
  test_plot <- plot_clusters(
    insert_colors(df_full_join(test_cluster_lists),5))
  expect_doppelganger("plot_API_all", test_plot)
})

test_that("plot_API() works", {
  test_plot <- plot_API(
    test_radii, test_centroids, 5,
    progbar = F, repulse = F
  )
  expect_doppelganger("plot_API_all", test_plot)
})

test_that("plot_API(rad_decrease = 0.8) works", {
  test_plot <- plot_API(
    test_radii, test_centroids, 5,
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
    test_radii, test_centroids, 5,
    progbar = F, repulse = F,
    try_place = T
  )
  expect_doppelganger(
    "plot_API_try_place_TRUE",
    test_plot
  )
})

# need to test more edge cases
