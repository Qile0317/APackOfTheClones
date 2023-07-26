source("testdata/SCIP.R")

test_that("run_apotc_warn_str works", {
	expect_identical(
		run_apotc_warn_str(test_integrated_pbmc, "_umap", NULL, 1),
		"No _umap reduction found on the seurat object, ensure the the reduction has been computed. Otherwise, did you mean umap ?"
	)
})
