test_that("defunct functions are defunct", {
	options(lifecycle_verbosity = "quiet")
	lifecycle::expect_defunct(integrate_tcr())
	lifecycle::expect_defunct(clonal_expansion_plot())
	lifecycle::expect_defunct(count_clone_sizes())
})
