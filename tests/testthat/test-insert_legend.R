# TODO
# sourcedata("v0", "cluster_lists")

# # assume correct
# plt <- plot_clusters(
# 	insert_colors(df_full_join(list(c1, c2, c3)), 3)
# )

# test_that("insert_legend() works", {
# 	expect_doppelganger(
# 		"plt_default_legend",
# 		insert_legend(plt, 1, 0)
# 	)

# 	expect_doppelganger(
# 		"plt_custom_legend",
# 		insert_legend(
# 			plt, 1, 0, c(0.5,2,17), "bottom_right", buffer = 2, #color = "black",
# 			legend_label = "custom label", legend_textsize = 4
# 		)
# 	)
# })

quietly_test_that("removeLegend works", {
    data("combined_pbmc")
    expect_equal(
        secretly_init_name(vizAPOTC(combined_pbmc, add_size_legend = FALSE)$layers),
        removeLegend(vizAPOTC(combined_pbmc))$layers
    )

    # TODO test if add layer
})

# quietly_test_that("overlayLegend works", {
#     data("combined_pbmc")
#     expect_equal(
#         vizAPOTC(combined_pbmc)$layers,
#         unname(overlayLegend(vizAPOTC(combined_pbmc))$layers)
#     )

#     # TODO test if add layer
# })
