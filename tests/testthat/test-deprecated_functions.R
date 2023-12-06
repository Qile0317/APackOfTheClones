source("testdata/SCIP.R")

test_that("integrate_tcr() works", {
	capture_output(
		integration_attempt <- suppressMessages(
			dev_integrate_tcr(test_pbmc, test_tcr, TRUE, FALSE)
		)
	)
	expect_identical(
		integration_attempt,
		test_integrated_pbmc
	)
})

test_that("integrate_tcr() works with verbose = FALSE", {
	expect_identical(
		dev_integrate_tcr(test_pbmc, test_tcr, verbose = FALSE, FALSE),
		test_integrated_pbmc
	)
})

test_that("percent_na() works", {
	expect_equal(percent_na(test_integrated_pbmc), 0)
	temp_pbmc <- test_integrated_pbmc
	for (i in 1:3) {temp_pbmc@meta.data[["barcode"]][[i]] <- NA}
	expect_equal(round(percent_na(temp_pbmc)), 4)
})

test_that("The main API works, where retain_axis_scales = F", {
	expect_doppelganger(
		"main_plot",
		clonal_expansion_plot(
			test_pbmc, test_tcr, verbose = FALSE, clone_scale_factor = 0.5,
			retain_axis_scales = FALSE, add_size_legend = FALSE
		)
	)
})

test_that("API can take an integrated object, retain_axis_scales = F", {
	expect_doppelganger(
		"main_plot",
		clonal_expansion_plot(
			test_integrated_pbmc, verbose = FALSE, clone_scale_factor = 0.5,
			retain_axis_scales = FALSE, add_size_legend = FALSE
		)
	)
})

test_that("API can use tsne and pca reductions, retain_axis_scales = F", {
	expect_doppelganger(
		"main_pca_plot",
		clonal_expansion_plot(
			test_integrated_pbmc, reduction = 'pca', verbose = FALSE,
			clone_scale_factor = 0.5, retain_axis_scales = FALSE,
			add_size_legend = FALSE
		)
	)
	expect_doppelganger(
		"main_tsne_plot",
		clonal_expansion_plot(
			test_integrated_pbmc, reduction = 'tsne', verbose = FALSE,
			clone_scale_factor = 0.5, retain_axis_scales = FALSE,
			add_size_legend = FALSE
		)
	)
})

test_that("retain_scale() works", {
	expect_doppelganger(
		"scale_retained_main_plot",
		suppressMessages(invisible(
			clonal_expansion_plot(
				test_integrated_pbmc, verbose = FALSE, clone_scale_factor = 0.5,
				retain_axis_scales = TRUE, add_size_legend = FALSE
			)
		))
	)
})

test_that("repulsion and try_place works", {
	expect_doppelganger(
		"try_place_plot",
		clonal_expansion_plot(
			test_integrated_pbmc, verbose = FALSE, repulse = FALSE,
			clone_scale_factor = 0.5, retain_axis_scales = FALSE,
			try_place = TRUE,
			add_size_legend = TRUE, legend_buffer = 1.5, legend_spacing = 0.25,
			legend_position = "bottom_right", legend_sizes = c(1, 5)
		)
	)
	expect_doppelganger(
		"repulsed_tight_plot",
		clonal_expansion_plot(
			test_integrated_pbmc, verbose = FALSE, repulse = TRUE,
			clone_scale_factor = 0.5, retain_axis_scales = FALSE,
			try_place = TRUE,
			add_size_legend = TRUE, legend_buffer = 1.5, legend_spacing = 0.25
		)
	)
})

# incredibly weird, try_place = TRUE works here but not on plot_API... WHAT?!?!
