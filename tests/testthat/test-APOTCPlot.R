# default args tested in test-noscript-defaultApotcPlot.R

test_that("add_default_theme works", {
	plt <- ggplot2::ggplot(data.frame("x" = c(1, 2), "y" = c(3, 4))) +
		ggplot2::geom_point(apotc_aes_string(x = "x", y = "y"))
	expect_doppelganger("add_default_theme_plot", add_default_theme(plt, "pca"))
})

# TODO more tests