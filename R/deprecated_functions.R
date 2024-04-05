md_deprecation_docstring <- function() {
	paste(
		'***ALL v0.1.x functions are deprecated***, and the workflow has been ',
		'completely revamped - now depending on the scRepertoire v2 package - ',
		'which allows for the processing of multi-sampled single cell data. ',
		'Please read the vignettes with `browseVignettes("APackOfTheClones")`',
		'or visit https://qile0317.github.io/APackOfTheClones/',
		sep = ""
	)
}

deprecation_docstring <- function() {
	gsub("*", "", md_deprecation_docstring(), fixed = TRUE)
}

#' DEFUNCT: Integrate a single TCR library into Seurat object metadata
#'
#' @description
#' `r lifecycle::badge("defunct")`
#'
#' `r md_deprecation_docstring()`
#'
#' @param ... arbitrary arguments
#'
#' @keywords internal
#'
#' @return error message
#' @export
#'
integrate_tcr <- function(...) {
	lifecycle::deprecate_stop(
		when = "1.0.0",
		what = I("integrate_tcr"),
		with = I("scRepertoire::combineExpression")
	)
}

#' @title
#' DEFUNCT: count the number of clonotype sizes per cell cluster in a seurat
#' object integrated with a TCR library
#'
#' @description
#' `r lifecycle::badge("defunct")`
#'
#' `r md_deprecation_docstring()`
#'
#' @param ... arbitrary arguments
#'
#' @keywords internal
#'
#' @return error message
#' @export
#'
count_clone_sizes <- function(...) {
	lifecycle::deprecate_stop(
		when = "1.0.0",
		what = I("`count_clone_sizes`"),
		with = I("countCloneSizes")
	)
}

#' @title
#' DEFUNCT: Visualize T cell clonal expansion with a ball-packing plot.
#'
#' @description
#' `r lifecycle::badge("defunct")`
#'
#' `r md_deprecation_docstring()`
#'
#' @param ... arbitrary arguments
#'
#' @keywords internal
#'
#' @return error message
#' @export
#'
clonal_expansion_plot <- function(...) {
	lifecycle::deprecate_stop(
		when = "1.0.0",
		what = I(
			"visualizing clonal expansion with `clonal_expansion_plot`"
		),
		with = I("vizAPOTC")
	)
}
