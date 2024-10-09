test_that("subsetting vizAPOTC works", {
	data("combined_pbmc")

	expect_doppelganger(
		"sample19",
		vizAPOTC(
			combined_pbmc,
			orig.ident = c("P19B", "P19L"),
			verbose = FALSE
		)
	)

	set.seed(829)
	expect_doppelganger(
		"void_retscale_tryplace_no_order_brl_c4:5",
		vizAPOTC(
			combined_pbmc,
			seurat_clusters = 4:5,
			retain_axis_scales = TRUE,
			try_place = TRUE,
			order_clones = FALSE,
			use_default_theme = FALSE,
			legend_position = "bottom right",
			add_legend_background = FALSE,
			verbose = FALSE
		)
	)
})
