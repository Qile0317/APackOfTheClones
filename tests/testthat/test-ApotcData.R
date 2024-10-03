data("combined_pbmc")
non_subset_apotc_data <- ApotcData(
	seurat_obj = combined_pbmc,
	alt_ident = NULL,
	metadata_filter_condition = "",
	clonecall = "CTstrict",
	reduction_base = "umap",
	clone_scale_factor = 0.300051,
	rad_scale_factor = 0.95
)

test_that("have_default_idents work", {
	expect_true(have_default_idents(non_subset_apotc_data))
})

expect_common_apotc_els_equal <- function(test_apotc_data) {

	expect_equal(
		test_apotc_data@reduction_base, non_subset_apotc_data@reduction_base
	
	)
	expect_equal(test_apotc_data@clonecall, non_subset_apotc_data@clonecall)
	
	expect_equal(get_num_clusters(test_apotc_data), get_num_clusters(non_subset_apotc_data))

	expect_equal(
		test_apotc_data@clone_scale_factor, non_subset_apotc_data@clone_scale_factor
	)
	expect_equal(
		test_apotc_data@rad_scale_factor, non_subset_apotc_data@rad_scale_factor
	)
	expect_equal(
		test_apotc_data@cluster_colors, non_subset_apotc_data@cluster_colors
	)
	expect_equal(test_apotc_data@labels, non_subset_apotc_data@labels)

}

test_that("The default case ApotcData constructor works", {
	test_apotc_data <- non_subset_apotc_data

	expected_centroids <- getdata("combined_pbmc", "all_cluster_centroids")

	expect_s4_class(test_apotc_data, "ApotcData")

	expect_equal(test_apotc_data@reduction_base, "umap")
	expect_equal(test_apotc_data@clonecall, "CTstrict")
	expect_equal(test_apotc_data@metadata_filter_string, "")

	expect_true(is_empty(test_apotc_data@clusters))

	expect_equal(test_apotc_data@centroids, expected_centroids)
	expect_equal(test_apotc_data@label_coords, expected_centroids)

	expect_equal(
		test_apotc_data@clone_sizes,
		getdata("get_clone_sizes", "raw_strict_clone_sizes")
	)

	expect_equal(get_num_clusters(test_apotc_data), 17)

	expect_equal(test_apotc_data@clone_scale_factor, 0.300051)
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

})

test_that("changing active ident affects clustering", {

	test_apotc <- ApotcData(
		seurat_obj = combined_pbmc,
		alt_ident = "seurat_clusters",
		metadata_filter_condition = "",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.300051,
		rad_scale_factor = 0.95
	)

	expect_common_apotc_els_equal(test_apotc)


})

test_that("The subset case ApotcData constructor works", {

	test_apotc_data <- ApotcData(
		seurat_obj = combined_pbmc,
		alt_ident = NULL,
		metadata_filter_condition = "seurat_clusters != 1",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.300051,
		rad_scale_factor = 0.95
	)

	expect_s4_class(test_apotc_data, "ApotcData")
	expect_true(is_empty(test_apotc_data@clusters))

	expect_common_apotc_els_equal(test_apotc_data)

	expect_equal(test_apotc_data@metadata_filter_string, "seurat_clusters != 1")

	expected_centroids <- getdata("combined_pbmc", "all_cluster_centroids")
	expected_centroids[[1]] <- list()
	expect_equal(test_apotc_data@centroids, expected_centroids)
	expect_equal(test_apotc_data@label_coords, expected_centroids)

	expected_raw_clone_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")
	expected_raw_clone_sizes[[1]] <- create_empty_table()
	expect_equal(
		test_apotc_data@clone_sizes, expected_raw_clone_sizes
	)

})

test_that("ApotcData subsetting works for 1 seurat cluster", {

	test_apotc_data <- ApotcData(
		seurat_obj = combined_pbmc,
		alt_ident = NULL,
		metadata_filter_condition = "seurat_clusters == 4",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.300051,
		rad_scale_factor = 0.95
	)

	expect_s4_class(test_apotc_data, "ApotcData")
	expect_true(is_empty(test_apotc_data@clusters))

	expect_common_apotc_els_equal(test_apotc_data)

	expect_equal(test_apotc_data@metadata_filter_string, "seurat_clusters == 4")

	expected_clone_sizes <- init_list(17, create_empty_table())
	expected_clone_sizes[[4]] <- getdata(
		"get_clone_sizes", "raw_strict_clone_sizes"
	)[[4]]
	expect_equal(
		object = get_raw_clone_sizes(test_apotc_data),
		expected = expected_clone_sizes
	)

})

test_that("circlepackClones works for the default case", {
	test_apotc_data <- circlepackClones(
		non_subset_apotc_data,
		ORDER = TRUE,
		try_place = FALSE,
		verbose = FALSE
	)

	expect_s4_class(test_apotc_data, "ApotcData")
	expect_common_apotc_els_equal(test_apotc_data)

	expect_equal(test_apotc_data@metadata_filter_string, "")

	expected_centroids <- getdata("combined_pbmc", "all_cluster_centroids")
	expected_centroids[[17]] <- list()

	expect_equal(test_apotc_data@centroids, expected_centroids)
	expect_equal(test_apotc_data@label_coords, expected_centroids)

	expect_equal(test_apotc_data@clone_sizes, non_subset_apotc_data@clone_sizes)

})

test_that("circlepackClones packs right for the default case", {
	skip_on_cran()

	test_apotc_data <- circlepackClones(
		non_subset_apotc_data,
		ORDER = TRUE,
		try_place = FALSE,
		verbose = FALSE
	)

	expect_equal(
		test_apotc_data@clusters, getdata("combined_pbmc", "expected_clusterlists")
	)
})

test_that("circlepackClones works for the subset case", {

	test_apotc_data <- ApotcData(
		seurat_obj = combined_pbmc,
		alt_ident = NULL,
		metadata_filter_condition = "seurat_clusters != 1",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.300051,
		rad_scale_factor = 0.95
	)

	test_apotc_data <- circlepackClones(
		test_apotc_data,
		ORDER = TRUE,
		try_place = FALSE,
		verbose = FALSE
	)

	expect_s4_class(test_apotc_data, "ApotcData")
	expect_common_apotc_els_equal(test_apotc_data)

	expect_equal(test_apotc_data@metadata_filter_string, "seurat_clusters != 1")

	expected_centroids <- getdata("combined_pbmc", "all_cluster_centroids")
	expected_centroids[[1]] <- list()
	expected_centroids[[17]] <- list()

	expect_equal(test_apotc_data@centroids, expected_centroids)
	expect_equal(test_apotc_data@label_coords, expected_centroids)

	expected_raw_clone_sizes <- getdata(
		"get_clone_sizes", "raw_strict_clone_sizes"
	)
	expected_raw_clone_sizes[[1]] <- create_empty_table()
	expect_equal(
		test_apotc_data@clone_sizes, expected_raw_clone_sizes
	)
	
})

test_that("circlepackClones packs right for the subset case", {
	skip_on_cran()
	
	test_apotc_data <- ApotcData(
		seurat_obj = combined_pbmc,
		alt_ident = NULL,
		metadata_filter_condition = "seurat_clusters != 1",
		clonecall = "CTstrict",
		reduction_base = "umap",
		clone_scale_factor = 0.300051,
		rad_scale_factor = 0.95
	)

	test_apotc_data <- circlepackClones(
		test_apotc_data,
		ORDER = TRUE,
		try_place = FALSE,
		verbose = FALSE
	)

	expected_clusterlists <- getdata("combined_pbmc", "expected_clusterlists")
	expected_clusterlists[[1]] <- list()
	expect_equal(test_apotc_data@clusters, expected_clusterlists)
})



test_that("match_clonotypes_to_sizes works", {

	a <- getLastApotcData(RunAPOTC(combined_pbmc, verbose = FALSE))

	expect_equal(
		match_clonotypes_to_sizes(a, get_unique_clonotypes(a)),
		get_aggregated_clone_sizes(a)
	)
})

test_that("match_index works", {

	expect_identical(match_index(non_subset_apotc_data, 1:17), 1:17)

	expect_identical(
		match_index(non_subset_apotc_data, index = gen_labels(17)), 1:17
	)

})



