source("testdata/cluster_lists.R")

test_that("df_full_join() works", {
  expect_equal(df_full_join(
    list(c1,c1_shifted_by_4_5, c1_shifted_to_9_0, c2,c3)),
    readRDS("testdata/cluster_df.rds"))
})
