source("testdata/cluster_lists.R")
suppressPackageStartupMessages(library(dplyr))

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
  expect_doppelganger(
    test_plot,
    "plot_clusters_c1_c2_c3")
})
