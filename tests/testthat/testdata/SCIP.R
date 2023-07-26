data("mini_clonotype_data","mini_seurat_obj")

test_pbmc <- mini_seurat_obj
test_tcr <- mini_clonotype_data
test_integrated_pbmc <- readRDS("testdata/mini_integrated_seurat_obj.rds")

test_apotc_obj_with_raw_clone_sizes <- new("apotc", clusters = list(NULL, NULL), num_clusters = 2, centroids = list(
	NULL, NULL), clone_sizes = list(list(structure(c(clonotype10 = 1L,
													 clonotype12 = 2L, clonotype15 = 1L, clonotype17 = 1L, clonotype18 = 1L,
													 clonotype19 = 1L, clonotype2 = 2L, clonotype20 = 2L, clonotype22 = 3L,
													 clonotype23 = 1L, clonotype24 = 4L, clonotype25 = 2L, clonotype26 = 2L,
													 clonotype27 = 4L, clonotype29 = 2L, clonotype3 = 2L, clonotype30 = 2L,
													 clonotype32 = 1L, clonotype33 = 2L, clonotype34 = 1L, clonotype35 = 2L,
													 clonotype37 = 3L, clonotype38 = 1L, clonotype39 = 1L, clonotype4 = 2L,
													 clonotype7 = 2L, clonotype8 = 1L, clonotype9 = 1L), dim = 28L, dimnames = list(
													 	x = c("clonotype10", "clonotype12", "clonotype15", "clonotype17",
													 		  "clonotype18", "clonotype19", "clonotype2", "clonotype20",
													 		  "clonotype22", "clonotype23", "clonotype24", "clonotype25",
													 		  "clonotype26", "clonotype27", "clonotype29", "clonotype3",
													 		  "clonotype30", "clonotype32", "clonotype33", "clonotype34",
													 		  "clonotype35", "clonotype37", "clonotype38", "clonotype39",
													 		  "clonotype4", "clonotype7", "clonotype8", "clonotype9")), class = "table")),
									list(structure(c(clonotype14 = 2L, clonotype16 = 1L, clonotype17 = 2L,
													 clonotype21 = 1L, clonotype24 = 2L, clonotype25 = 1L, clonotype27 = 1L,
													 clonotype3 = 3L, clonotype32 = 1L, clonotype33 = 1L, clonotype35 = 2L,
													 clonotype36 = 2L, clonotype37 = 1L, clonotype38 = 2L, clonotype39 = 1L,
													 clonotype5 = 2L, clonotype7 = 1L, clonotype8 = 1L, clonotype9 = 3L
									), dim = 19L, dimnames = list(x = c("clonotype14", "clonotype16",
																		"clonotype17", "clonotype21", "clonotype24", "clonotype25",
																		"clonotype27", "clonotype3", "clonotype32", "clonotype33",
																		"clonotype35", "clonotype36", "clonotype37", "clonotype38",
																		"clonotype39", "clonotype5", "clonotype7", "clonotype8",
																		"clonotype9")), class = "table"))), clone_scale_factor = 1,
	rad_scale_factor = 1, cluster_colors = c("#F8766D", "#00BFC4"
	), reduction_base = "umap", labels = character(0), label_coords = list())
