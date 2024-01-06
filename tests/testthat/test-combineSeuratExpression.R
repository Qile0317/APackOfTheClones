test_that("combineSeuratExpression output is compatible with scRepertoire", {

	# scRepertoire data
	data("scRep_example", "mini_contig_list", package = "scRepertoire")

	# named test
	mini_combined <- scRepertoire::combineTCR(
		input.data = scRepertoire::mini_contig_list,
		samples = c("P17B", "P17L", "P18B", "P18L",
					"P19B","P19L", "P20B", "P20L")
	)
	apotc_pbmc <- combineSeuratExpression(
		mini_combined, scRepertoire::scRep_example
	)
	screp_pbmc <- scRepertoire::combineExpression(
		mini_combined, scRepertoire::scRep_example, cloneCall = "strict",
		chain = "both", group.by = NULL, proportion = FALSE, filterNA = FALSE,
		cloneSize = c(
			Rare = 1, Small = 5, Medium = 50, Large = 100, Hyperexpanded = 1000
		),
		addLabel = FALSE
	)

	apotc_pbmc@commands[[length(apotc_pbmc@commands)]] <- NULL
	screp_pbmc@commands[[length(screp_pbmc@commands)]] <- NULL

	expect_equal(apotc_pbmc, screp_pbmc)

	# TODO unnamed test! (was a problem in scRepertoire)
})

# TODO try this on the v0 objects
