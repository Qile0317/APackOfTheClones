source("testdata/cluster_lists.R")

test_that("df_full_join() works", {
    test_df <- readRDS("testdata/df_full_join_test_df.rds")

    expect_equal(
        df_full_join(
            list(c1, c1_shifted_by_4_5, c1_shifted_to_9_0, c2,c3)
        ),
        test_df
    )

    expect_equal(
        df_full_join(list(c1)),
        test_df[1:9, ]
    )

    expect_equal(df_full_join(list(c1, NA)), test_df[1:9, ])
    expect_equal(df_full_join(list(c1, list())), test_df[1:9, ])

    test_df_rows_1_till_9_with_label_1 <- test_df[1:9,]
    test_df_rows_1_till_9_with_label_1[[1]] <- rep("cluster 1", 9)

    expect_equal(
        df_full_join(list(NA, c1)),
        test_df_rows_1_till_9_with_label_1
    )

    expect_equal(
        test_df_rows_1_till_9_with_label_1,
        df_full_join(list(list(), c1, NA))
    )

    expect_equal(
        df_full_join(list(NA, c1, NA, c2, NA)),
        df_full_join(list(list(), c1, list(), c2, list()))
    )
})
# NA testcases shouldnt be needed anymore (?)

test_that("plot_clusters() works", {
    test_plot <- plot_clusters(
        insert_colors(df_full_join(list(c1,c2,c3)), 3))

    expect_doppelganger("plot_c1_c2_c3", test_plot)

    test_plot <- plot_clusters(
        insert_colors(df_full_join(test_cluster_lists), 5))

    expect_doppelganger("plot_API_all", test_plot)
})

test_that("plot_API() works", {
    test_plot <- plot_API(
        test_radii, test_centroids, 5, try_place = FALSE,
        progbar = FALSE, repulse = FALSE
    )
    expect_doppelganger("plot_API_all", test_plot)
})

test_that("plot_API(rad_decrease = 0.2) works", {
    test_plot <- plot_API(
        test_radii, test_centroids, 5, try_place = FALSE,
        progbar = FALSE, repulse = FALSE, rad_decrease = 0.8
    )
    expect_doppelganger(
        "plot_API_all_decreased_0.2",
        test_plot
    )
})

# the following test produces inconsistent plots locally. however, when used in
# clonal_expansion_plot it works as intended. All its member functions seem to
# also work as intended though it hasnt been examined in detail
# in a future patch the cause of this will be attempted to be fixed.

# i have a hypothesis - this might be failing in the r cmd check
# environment because there sometimes might be multiple best ways
# to place "as close as possible"? As evidenced by previous tests too
#
# it fails when running devtools::test() on opening RStudio, but not a bit after
# when more packages are stuff are loading in. It might be some hidden dumb dep
# problem that alters the behaviour of a function? No idea. Or maybe floats
# being slightly differenty valued?
#
# the changed snapshot shows that bottom two circles being differently placed
# with all else equal

#test_that("plot_API(try_place = TRUE) works", {
#  test_plot <- plot_API(
#    test_radii, test_centroids, 5, ORDER = TRUE, try_place = TRUE,
#    progbar = FALSE, repulse = FALSE
#  )
#  expect_doppelganger(
#    "plot_API_try_place",
#    test_plot
#  )
#})
