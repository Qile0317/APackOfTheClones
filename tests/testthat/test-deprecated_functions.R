sourcedata("deprecated_functions", "SCIP")
test_integrated_pbmc <- getdata(
    "deprecated_functions", "mini_integrated_seurat_obj"
)

test_that("integrate_tcr() works and is deprecated", {

	options(lifecycle_verbosity = "quiet")
	lifecycle::expect_deprecated(
		integrate_tcr(test_pbmc, test_tcr, verbose = FALSE)
	)

	capture_output(
		integration_attempt <- suppressMessages(
			dev_integrate_tcr(test_pbmc, test_tcr, verbose = TRUE)
		)
	)
	expect_identical(
		integration_attempt, test_integrated_pbmc
	)

	expect_identical(
		dev_integrate_tcr(test_pbmc, test_tcr, verbose = FALSE),
		test_integrated_pbmc
	)
})

test_that("percent_na() works", {
	expect_equal(percent_na(test_integrated_pbmc), 0)
	temp_pbmc <- test_integrated_pbmc
	for (i in 1:3) temp_pbmc@meta.data[["barcode"]][[i]] <- NA
	expect_equal(round(percent_na(temp_pbmc)), 4)
})

test_that("count_clone_sizes works and is deprecated", {
	options(lifecycle_verbosity = "quiet")
	lifecycle::expect_deprecated(
		count_clone_sizes(test_integrated_pbmc)
	)
	trial <- count_clone_sizes(test_integrated_pbmc)

	expect_equal(length(trial), 2)

	untable <- function(a) as.numeric(unname(a))
	expect_equal(untable(trial[[1]]), c(12, 12, 2, 2))
	expect_equal(untable(trial[[2]]), c(10, 7, 2))
})

test_that("clonal_expansion_plot is defunct", {
	options(lifecycle_verbosity = "quiet")
	lifecycle::expect_defunct(clonal_expansion_plot())
})
