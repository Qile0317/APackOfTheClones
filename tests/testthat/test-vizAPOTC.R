test_that("vizAPOTC works", {
    data("combined_pbmc")
	expect_doppelganger(
		.defaultApotcPlot,
		vizAPOTC(combined_pbmc,  verbose = FALSE)
	)
})

test_that("subsetting vizAPOTC works", {
	data("combined_pbmc")

	expect_doppelganger(
		.sample19ApotcPlot,
		vizAPOTC(
			combined_pbmc,
			orig.ident = c("P19B", "P19L"),
			verbose = FALSE
		)
	)

	set.seed(829)
	expect_doppelganger(
		"ret_scale_try_place_scrambled_5",
		vizAPOTC(
			combined_pbmc,
			seurat_clusters = 1:5,
			retain_axis_scales = TRUE,
			try_place = TRUE,
			order_clones = FALSE,
			verbose = FALSE
		)
	)
})