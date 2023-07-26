source("testdata/SCIP.R")

test_that("get_legend_coordinate works", {
	plt <- plot_clusters(df_full_join(list(c1,c2,c3)))

	expect_equal(get_legend_coordinate(plt, c(1,2), 1.5), c(1, 2))
	expect_equal(get_legend_coordinate(plt, "top_left", 1.5),
})

test_that("insert_legend() works", {
	plt <- plot_clusters(df_full_join(list(c1,c2,c3)))
	expect_doppelganger(
		"plt_default_legend",
		insert_legend(plt, 1, 0)
	)

	expect_doppelganger(
		"plt_custom_legend",
		insert_legend(
			plt, 1, 0, c(0.5,2,17), "bottom_right", buffer = 2, #color = "black",
			legend_label = "custom label", legend_textsize = 4
		)
	)
})
