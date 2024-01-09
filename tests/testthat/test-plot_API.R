sourcedata("v0", "cluster_lists")

test_that("df_full_join() works", {
    test_df <- getdata("plot_API", "df_full_join_test_df")

    expect_equal(
        df_full_join(
            list(c1, c1_shifted_by_4_5, c1_shifted_to_9_0, c2, c3)
        ),
        test_df
    )

    expect_equal(
        df_full_join(list(c1)),
        test_df[1:9, ]
    )

    expect_equal(df_full_join(list(c1, list())), test_df[1:9, ])

    test_df_rows_1_till_9_with_label_1 <- test_df[1:9,]
    test_df_rows_1_till_9_with_label_1[[1]] <- rep("cluster 1", 9)

    expect_equal(
        df_full_join(list(list(), c1)),
        test_df_rows_1_till_9_with_label_1
    )

    expect_equal(
        test_df_rows_1_till_9_with_label_1,
        df_full_join(list(list(), c1, list()))
    )

})

test_that("plot_clusters() works", {
    test_plot <- plot_clusters(
        insert_colors(df_full_join(list(c1, c2, c3)), 3))

    expect_doppelganger("plot_c1_c2_c3", test_plot)

    test_plot <- plot_clusters(
        insert_colors(df_full_join(test_cluster_lists), 5))

    expect_doppelganger("plot_API_all", test_plot)
})

# TODO more tests
