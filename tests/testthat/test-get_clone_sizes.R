# TODO need to redo everything here, based on filtering confition
# NOT COMPREHENSIVE

data("combined_pbmc")
expected_clone_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")

test_that("count_raw_clone_sizes works", {
	expect_identical(
		count_raw_clone_sizes(get(data("combined_pbmc")), 17, "CTstrict"),
		expected_clone_sizes
	)
})

test_that("countCloneSizes works", {

	expect_identical(countCloneSizes(combined_pbmc), expected_clone_sizes)

	expect_identical(
		object = countCloneSizes(combined_pbmc, seurat_clusters = 1),
		expected = append(
			list(expected_clone_sizes[[1]]),
			init_list(length(expected_clone_sizes) - 1, create_empty_table())
		)
	)

})

test_that("mergeCloneSizes works", {

	expect_identical(
		mergeCloneSizes(countCloneSizes(combined_pbmc)),
		mergeCloneSizes(expected_clone_sizes)
	)

	test_merged_cluster1 <- mergeCloneSizes(
		countCloneSizes(combined_pbmc, seurat_clusters = 1),
		sort_decreasing = TRUE
	)

	expect_contains(test_merged_cluster1, expected_clone_sizes[[1]])

})

quietly_test_that("get_top_clonotypes works", {
	expect_contains(
		get_top_clonotypes(expected_clone_sizes, Inf),
		get_unique_clonotypes(getLastApotcData(RunAPOTC(get(data("combined_pbmc")))))
	)
	#TODO - not remotely a good test
})

# TODO other functions
