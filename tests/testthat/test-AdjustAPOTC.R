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
		suppressMessages(APOTCPlot(AdjustAPOTC(combined_pbmc, verbose = FALSE)))
	)

  # change clone and rad scale, nudge first clusters 1 & 2,
  # move clusters 3 & 4, recolor 5th and 6th clusters
  combined_pbmc <- AdjustAPOTC(combined_pbmc, run_id = "r2", verbose = FALSE,
    new_rad_scale_factor = 0.9,
    new_clone_scale_factor = 0.35,
    nudge_cluster = c(1, 2),
    nudge_vector = c(1, 1),
    relocate_cluster = c(3, 4),
    relocation_coord = list(c(-0.5, -2.5), c(3, -4.5)),
    recolor_cluster = c(5, 6),
    new_color = c("green", "cyan")
  )

  expect_doppelganger(
		"AdjustAPOTC_plot",
		suppressMessages(APOTCPlot(combined_pbmc, run_id = "r2"))
	)

  # sanity check - r1 isnt modified
  expect_doppelganger(
		.defaultApotcPlot,
		APOTCPlot(
      AdjustAPOTC(combined_pbmc, run_id = "r1", verbose = FALSE),
      run_id = "r1"
    )
	)

  # repulse and add labels
  combined_pbmc <- AdjustAPOTC(combined_pbmc, run_id = "r2", verbose = FALSE,
    repulse = TRUE
  )

  expect_doppelganger(
		"AdjustAPOTC_plot_repulsed_with_label",
    APOTCPlot(combined_pbmc, run_id = "r2", show_labels = TRUE, label_size = 7)
	)

})
