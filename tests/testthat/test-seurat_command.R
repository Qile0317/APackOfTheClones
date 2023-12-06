data("mini_seurat_obj", "mini_clonotype_data")

# doesnt rlly need testing as its just copied directly from Seurat

test_that("make_apotc_command works", {
	new_mini_seurat_obj <- RunAPOTC(
		mini_seurat_obj, mini_clonotype_data, clone_scale_factor = 0.1,
		verbose = FALSE
	)

	expect_equal(
		new_mini_seurat_obj@commands[["RunAPOTC"]]@params,
		list(
			tcr_df = "removed from command call to save memory",
			reduction_base = "umap", clone_scale_factor = 0.1,
			rad_scale_factor = 0.95, ORDER = TRUE, scramble = FALSE,
			try_place = FALSE, repulse = FALSE, repulsion_threshold = 1,
			repulsion_strength = 1, max_repulsion_iter = 10L, verbose = FALSE
		)
	)

	expect_identical(
		new_mini_seurat_obj@commands[["RunAPOTC"]]@name, "RunAPOTC"
	)
	expect_identical(
		new_mini_seurat_obj@commands[["RunAPOTC"]]@assay.used, "RNA"
	)
	expect_identical(
		new_mini_seurat_obj@commands[["RunAPOTC"]]@call.string,
		c(
			paste(
				"RunAPOTC(mini_seurat_obj, mini_clonotype_data,",
				 "clone_scale_factor = 0.1, "
			),
			"    verbose = FALSE)"
		)
	)
})

test_that("get_cmd works", {
	new_mini_seurat_obj <- RunAPOTC(
		mini_seurat_obj, mini_clonotype_data, verbose = FALSE
	)
	expect_identical(get_cmd(new_mini_seurat_obj, "ORDER"), TRUE)
})
