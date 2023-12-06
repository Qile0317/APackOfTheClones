test_that("combineSeuratExpression output is compatible with scRepertoire", {
	data("screp_example") # from scRepertoire
	data("combined_contigs")

	apotc_pbmc <- combineSeuratExpression(combined_contigs, screp_example)
	screp_pbmc <- combineExpression(combined_contigs, screp_example)

	screp_pbmc_meta_modified <- screp_pbmc@meta.data
	screp_pbmc_meta_modified[["barcode"]] <- NULL

	# compatible doesnt mean identical. In apotc:
	# - extra barcodes are removed from metadata
	# - there is a seurat command slot

	expect_equal(apotc_pbmc@assays, screp_pbmc@assays)
	expect_equal(apotc_pbmc@meta.data, screp_pbmc_meta_modified)
	expect_identical(apotc_pbmc@active.assay, screp_pbmc@active.assay)
	expect_identical(apotc_pbmc@active.ident, screp_pbmc@active.ident)
	expect_identical(apotc_pbmc@graphs, screp_pbmc@graphs)
	expect_identical(apotc_pbmc@neighbors, screp_pbmc@neighbors)
	expect_equal(apotc_pbmc@reductios, screp_pbmc@reductions)
	expect_identical(apotc_pbmc@images, screp_pbmc@images)
	expect_identical(apotc_pbmc@project.name, screp_pbmc@project.name)
	expect_identical(apotc_pbmc@misc, screp_pbmc@misc)
	expect_identical(apotc_pbmc@version, screp_pbmc@version)
	# commands should be seperately tested
	# tools is probably not too relevant
})
