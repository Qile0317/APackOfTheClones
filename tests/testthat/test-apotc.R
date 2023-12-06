source("testdata/SCIP.R")

#test_that("the apotc class works", {})

test_that("run_apotc_warn_str works", {
	expect_identical(
		run_apotc_warn_str(test_integrated_pbmc, "APOTC", 1, TRUE, FALSE),
		"No _umap reduction found on the seurat object, ensure the the reduction has been computed. Otherwise, did you mean umap ?"
	)
	expect_identical(
		run_apotc_warn_str(test_integrated_pbmc, "_umap", 1, TRUE, FALSE),
		"No _umap reduction found on the seurat object, ensure the the reduction has been computed. Otherwise, did you mean umap ?"
	)
	expect_identical(
		run_apotc_warn_str(test_integrated_pbmc, "umap", -0.01, TRUE, FALSE),
		"clone_scale_factor has to be a positive number"
	)
	expect_identical(
		run_apotc_warn_str(test_integrated_pbmc, "umap", 1, TRUE, TRUE),
		"ORDER and scramble are both TRUE, please set only one to TRUE"
	)
	expect_null(run_apotc_warn_str(
		test_integrated_pbmc, "umap", "auto", TRUE, FALSE
	))
})

test_that("RunAPOTC works", {

})
