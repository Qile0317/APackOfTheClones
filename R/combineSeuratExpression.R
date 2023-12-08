#' @title
#' Wrapper to add scRepertoire-compatible clonotype information to a seurat
#' object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This is a modified version of [scRepertoire::combineExpression], which adds
#' the immune receptor information to a Seurat object's meta data. This function
#' also calculates the combined clonotype frequencies.
#'
#' Importantly, before using this function, ensure the barcodes of the seurat
#' object match the barcodes in the output of the [scRepertoire::combineTCR]
#' or [scRepertoire::combineBCR] call. Check [AddSampleBarcode]
#' to change the prefixes of the Seurat object. If combining more than one
#' immune receptor type, barcodes with both receptors will be removed during the
#' combination process.
#'
#' @details
#' This function is a sort of "special case" of the original function
#' that combines and processes the clonotype data so that it can be used by
#' `APackOfTheClones` further downstream. The modified output seurat object
#' should be completely compatible to be used with other `scRepertoire`
#' functions.
#'
#' @param input.data The product of [scRepertoire::combineTCR],
#' [scRepertoire::combineBCR] or a list of both. (see scRepertoire's vignette)
#' @param sc.data The Seurat object to attach. Note that this will not work with
#' SingleCellExperiment (SCE) objects!
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa),
#' VDJC gene + CDR3 nucleotide (strict) or a custom variable in the data.
#' @param chain indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column label in the combined clones in which
#' clonotype frequency will be calculated. "none" will keep the list as is,
#' while NULL will merge all the clones into a single data frame.
#' @param proportion Whether to proportion (**TRUE**) or total frequency
#' (**FALSE**) of the clonotype based on the group.by variable. For
#' `APackOfTheClones`, this value is recommended (and defaults) to be set
#' *FALSE* for cleaner cluster visualization.
#' @param cloneSize The bins for the grouping based on proportion or frequency.
#' If proportion is **FALSE** and the *cloneSizes* are not set high enough based
#' on frequency, the upper limit of *cloneSizes* will be automatically amended.
#' @param filterNA Method to subset seurat object of barcodes without
#' clonotype information
#' @param addLabel This will add a label to the frequency header, allowing
#' the user to try multiple group.by variables or recalculate frequencies after
#' subsetting the data.
#'
#' @return
#' The modified seurat object with the contig information attached, and a seurat
#' command slot with the command call information
#'
#' @examples
#' #unfinished
#'
#' @export
#'
combineSeuratExpression <- function(
	input.data,
	sc.data,
	cloneCall = "strict",
	chain = "both",
	group.by = NULL,
	proportion = FALSE,
	filterNA = FALSE,
	cloneSize = c(
		Rare = 1, Small = 5, Medium = 50, Large = 100, Hyperexpanded = 1000
	),
	addLabel = FALSE
) {
	call_time <- Sys.time()

	# handle inputs
	if (!is_seurat_object(sc.data)) {
		stop("`sc.data` must be a Seurat object")
	}

	# run combineExpression
	combined_seurat <- scRepertoire::combineExpression(
		input.data = input.data,
		sc.data = sc.data,
		cloneCall = cloneCall,
		chain = chain,
		group.by = group.by,
		proportion = proportion,
		filterNA = filterNA,
		cloneSize = cloneSize,
		addLabel = addLabel
	)

	# overwrite the seurat command slot
	combined_seurat@commands[["combineExpression"]] <- NULL
	combined_seurat@commands[["combineSeuratExpression"]] <- make_apotc_command(
		call_time = call_time
	)

	combined_seurat
}

# probably should also have a function in APOTC to add a metadata slot for pure
# counts
