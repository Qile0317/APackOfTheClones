source("testdata/cluster_lists.R")
suppressPackageStartupMessages(library(dplyr))
library(vdiffr)
library(ggforce)

test_that("df_full_join() works", {
  test_df <- readRDS("testdata/cluster_df.rds")

  expect_equal(df_full_join(
    list(c1,c1_shifted_by_4_5, c1_shifted_to_9_0, c2,c3)),
    test_df)

  expect_equal(df_full_join(list(c1)),
               test_df[1:9, ])
})

test_that("plot_clusters() works", {
  cluster_df <- df_full_join(list(c1,c2,c3))
  test_plot <- plot_clusters(cluster_df)
  expect_doppelganger("plot_c1_c2_c3", test_plot)
})

test_that("plot_API() works", {
  test_plot <- plot_API(
    test_radii[-1], test_centroids[-1],
    progbar = F, repulse = F
  )
  expect_doppelganger("plot_API_all", test_plot)
})

test_that("plot_API(rad_decrease = 0.8) works", {
  test_plot <- plot_API(
    test_radii[-1], test_centroids[-1],
    progbar = F, repulse = F,
    rad_decrease = 0.8
  )
  expect_doppelganger(
    "plot_API_all_scaled_0.8",
    test_plot
  )
})

# need to test try_place
# need to test more edge cases
