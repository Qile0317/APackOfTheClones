test_that("AdjustAPOTC works", {
  combined_pbmc <- RunAPOTC(
    get(data("combined_pbmc")), run_id = "r1", verbose = FALSE
  )

  combined_pbmc <- RunAPOTC(
    combined_pbmc, run_id = "r2", verbose = FALSE
  )

  # sanity check - Adjust nothing
  expect_doppelganger(
		.defaultApotcPlot,
		APOTCPlot(AdjustAPOTC(combined_pbmc, verbose = FALSE))
	)

  # TODO adjust 

})
