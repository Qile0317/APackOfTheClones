# TODO need to redo everything here, based on filtering confition
# NOT COMPREHENSIVE

test_that("count_raw_clone_sizes works", {
	expect_identical(
		count_raw_clone_sizes(combined_pbmc, 17, "CTstrict"),
		getdata("get_clone_sizes", "raw_strict_clone_sizes")
	)
})

# TODO other functions
