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
})