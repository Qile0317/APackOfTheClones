test_that("vizAPOTC works", {
    data("combined_pbmc")
	expect_doppelganger(
		.defaultApotcPlot,
		vizAPOTC(combined_pbmc,  verbose = FALSE)
	)
})
