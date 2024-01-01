test_that("The default case ApotcData constructor works", {

	test_apotc_data <- ApotcData(
		seurat_obj = get(data("combined_pbmc")),
		metadata_filter_condition = "",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.95,
		rad_scale_factor = 0.95
	)

	expected_centroids <- getdata("combined_pbmc", "all_cluster_centroids")

	expect_s4_class(test_apotc_data, "ApotcData")

	expect_equal(test_apotc_data@reduction_base, "umap")
	expect_equal(test_apotc_data@clonecall, "CTstrict")
	expect_equal(test_apotc_data@metadata_filter_string, "")

	expect_false(isnt_empty(test_apotc_data@clusters)) # TODO test this??

	expect_equal(test_apotc_data@centroids, expected_centroids)
	expect_equal(
		test_apotc_data@clone_sizes,
		getdata("get_clone_sizes", "raw_strict_clone_sizes")
	)

	expect_equal(test_apotc_data@num_clusters, 17)

	expect_equal(test_apotc_data@clone_scale_factor, 0.95)
	expect_equal(test_apotc_data@rad_scale_factor, 0.95)
	expect_equal(test_apotc_data@cluster_colors, c(
		"#F8766D", "#E7851E", "#D09400", "#B2A100", "#89AC00", "#45B500",
		"#00BC51", "#00C087", "#00C0B2", "#00BCD6", "#00B3F2", "#29A3FF",
		"#9C8DFF", "#D277FF", "#F166E8", "#FF61C7", "#FF689E"
	))
	expect_equal(test_apotc_data@labels, c(
		"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
		"C13", "C14", "C15", "C16", "C17"
	))
	expect_equal(test_apotc_data@label_coords, expected_centroids)
})

test_that("The subset case ApotcData constructor works", {

	test_apotc_data <- ApotcData(
		seurat_obj = get(data("combined_pbmc")),
		metadata_filter_condition = "seurat_clusters != 1",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.95,
		rad_scale_factor = 0.95
	)

	expect_s4_class(test_apotc_data, "ApotcData")

	expect_equal(test_apotc_data@reduction_base, "umap")
	expect_equal(test_apotc_data@clonecall, "CTstrict")
	expect_equal(test_apotc_data@metadata_filter_string, "seurat_clusters != 1")

	expect_false(isnt_empty(test_apotc_data@clusters)) # TODO test this??

	expected_centroids <- getdata("combined_pbmc", "all_cluster_centroids")
	expected_centroids[[1]] <- list()
	expect_equal(test_apotc_data@centroids, expected_centroids)

	expected_raw_strict_clone_sizes <- getdata(
		"get_clone_sizes", "raw_strict_clone_sizes"
	)
	expected_raw_strict_clone_sizes[[1]] <- create_empty_table()
	expect_equal(
		test_apotc_data@clone_sizes, expected_raw_strict_clone_sizes
	)

	expect_equal(test_apotc_data@num_clusters, 17)

	expect_equal(test_apotc_data@clone_scale_factor, 0.95)
	expect_equal(test_apotc_data@rad_scale_factor, 0.95)
	expect_equal(test_apotc_data@cluster_colors, c(
		"#F8766D", "#E7851E", "#D09400", "#B2A100", "#89AC00", "#45B500",
		"#00BC51", "#00C087", "#00C0B2", "#00BCD6", "#00B3F2", "#29A3FF",
		"#9C8DFF", "#D277FF", "#F166E8", "#FF61C7", "#FF689E"
	))
	expect_equal(test_apotc_data@labels, c(
		"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
		"C13", "C14", "C15", "C16", "C17"
	))

	expect_equal(test_apotc_data@label_coords, expected_centroids)
})


# TODO test the circle packing with diff args, the @clusters slot shgould be identical to other tests

# TODO test the repulsion API

# TODO test getters and setters
