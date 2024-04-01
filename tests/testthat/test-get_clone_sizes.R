# TODO need to redo everything here, based on filtering confition
# NOT COMPREHENSIVE

expected_clone_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")

test_that("count_raw_clone_sizes works", {
	expect_identical(
		count_raw_clone_sizes(get(data("combined_pbmc")), 17, "CTstrict"),
		expected_clone_sizes
	)
})

test_that("countCloneSizes works", {
	expect_identical(
		countCloneSizes(get(data("combined_pbmc"))),
		expected_clone_sizes
	)
})

quietly_test_that("get_top_clonotypes works", {
	expect_contains(
		get_top_clonotypes(expected_clone_sizes, Inf),
		get_unique_clonotypes(getLastApotcData(RunAPOTC(get(data("combined_pbmc")))))
	)
	#TODO - not remotely a good test
})

# TODO other functions
