
# NOT COMPREHENSIVE

data("combined_pbmc")
expected_clone_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")

test_that("count_raw_clone_sizes works", {
	expect_identical(
		count_raw_clone_sizes(combined_pbmc, as.character(1:17), "CTstrict"),
		expected_clone_sizes
	)
})

test_that("countCloneSizes works", {

	expect_identical(unname(countCloneSizes(combined_pbmc)), expected_clone_sizes)

	expect_identical(
		object = unname(countCloneSizes(combined_pbmc, seurat_clusters = 1)),
		expected = append(
			list(expected_clone_sizes[[1]]),
			init_empty_table_list(length(expected_clone_sizes) - 1)
		)
	)

	expect_contains(
        object = combined_pbmc %>%
			countCloneSizes(extra_filter = "clonalFrequency > 1L") %>%
			unname() %>%
			unlist() %>%
			names(),
        expected = names(getdata("ApotcClonalNetwork", "shared_clones"))
    )

	same_ident_pbmc <- combined_pbmc
	same_ident_pbmc@active.ident <- "Foo" %>%
		rep(length(combined_pbmc@active.ident)) %>%
		as.factor()

	test_obj <- countCloneSizes(same_ident_pbmc, sort_decreasing = TRUE)[[1]]
	expected_obj <- countCloneSizes(
		combined_pbmc, by_cluster = FALSE, sort_decreasing = TRUE
	)

	expect_mapequal(test_obj, expected_obj)
	expect_identical(test_obj[1], expected_obj[1]) # size 11
	expect_mapequal(test_obj[2:4], expected_obj[2:4]) # size 3
	expect_mapequal(test_obj[5:12], expected_obj[5:12]) # size 2

	expect_identical(
		unname(countCloneSizes(combined_pbmc, by_cluster = "seurat_clusters")),
		expected_clone_sizes
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

	expect_identical(test_merged_cluster1, expected_clone_sizes[[1]])

})

quietly_test_that("get_top_clonotypes works", {

	expect_setequal(
		get_top_clonotypes(expected_clone_sizes, Inf),
		get_unique_clonotypes(getLastApotcData(RunAPOTC(combined_pbmc)))
	)


})

test_that("aggregate_clone_sizes works", {

	expect_identical(
		object = aggregate_clone_sizes(list()), expected = numeric(0)
	)

	test_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")

	expect_identical(
		object = aggregate_clone_sizes(test_sizes[1]),
		expected = convert_table_to_named_numeric(test_sizes[[1]])
	)

	expect_identical(
		object = aggregate_clone_sizes(list(table(letters)), sort_decreasing = TRUE),
		convert_table_to_named_numeric(table(letters))
	)

})


