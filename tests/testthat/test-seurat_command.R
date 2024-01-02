test_that("everything in seurat_command works", {

	data("combined_pbmc")
	combined_pbmc <- RunAPOTC(combined_pbmc, verbose = FALSE)

	expect_identical(
		getlast(names(combined_pbmc@commands)),
		"RunAPOTC.umap;CTstrict;_;_"
	)

	test_cmd_obj <- getlast(combined_pbmc@commands)
	expect_s4_class(test_cmd_obj, "SeuratCommand")
	
	expect_equal(
		test_cmd_obj@params,
		list(reduction_base = "umap", clonecall = "CTstrict", clone_scale_factor = 0.300051, 
		rad_scale_factor = 0.95, order_clones = TRUE, scramble_clones = FALSE, 
		try_place = FALSE, repulse = TRUE, repulsion_threshold = 1, 
		repulsion_strength = 1, max_repulsion_iter = 20L, override = TRUE, 
		verbose = FALSE)
	)

	expect_identical(test_cmd_obj@name, "RunAPOTC")
	expect_identical(test_cmd_obj@assay.used, "RNA")
	expect_identical(
		test_cmd_obj@call.string, "RunAPOTC(combined_pbmc, verbose = FALSE)"
	)
})
