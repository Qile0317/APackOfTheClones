# all functions in this script are copied/modified from scRepertoire - with
# permission from the author

ce_warn_error_helper <- function() {
	#params <- get_parent_params(n = 2)
	#print(params) # unfinished, this is just a test
}

clonecall_to_colname <- function(x) {
	# need to make this a global var :P
	.clonecall_dictionary <- hash::hash(
		"ctnt" = "CTnt",
		"ctgene" = "CTgene",
		"ctaa" = "CTaa",
		"CTstrict" = "CTstrict",
		"nt" = "CTnt",
		"nucleotide" = "CTnt",
		"aa" = "CTaa",
		"amino" = "CTaa",
		"gene+nt" = "CTstrict",
		"strict" = "CTstrict"
	)

	if (x %in% c("CTnt", "CTgene", "CTaa", "CTstrict")) {return(x)}
	x <- tolower(x)
	if (!is.null(.clonecall_dictionary[[x]])) {
		return(.clonecall_dictionary[[x]])
	}
	stop(paste(
		"invalid input cloneCall, did you mean: '",
		closest_word(x, names(.clonecall_dictionary)),
		"'?",
		sep = ""
	))
}

# grab the metadata of a seurat object sc
grabMeta <- function(sc) {
	meta <- data.frame(sc[[]], slot(sc, "active.ident"))
	colnames(meta)[length(meta)] <- "ident"
	meta
}

off.the.chain <- function(dat, chain, cloneCall) {
	chain1 <- toupper(chain)
	if (chain1 %in% c("TRA", "TRG", "IGH")) {
		x <- 1
	} else if (chain1 %in% c("TRB", "TRD", "IGL")) {
		x <- 2
	} else {
		warning("It looks like ", chain, " does not match the available options for `chain = `")
	}
	dat[,cloneCall] <- stringr::str_split(dat[,cloneCall], "_", simplify = TRUE)[,x]
	return(dat)
}

filteringNA <- function(sc, cloneType) {
	meta <- grabMeta(sc)
	evalNA <- data.frame(meta[,"cloneType"])
	colnames(evalNA) <- "indicator"
	evalNA <- evalNA %>%
		dplyr::transmute(indicator = ifelse(is.na(indicator), 0, 1))
	rownames(evalNA) <- rownames(meta)

	col.name <- names(evalNA) %||% colnames(evalNA)
	sc[[col.name]] <- evalNA
	subset(sc, cloneType != 0)
}

add_clonetype_porportion <- function(Con.df, cloneTypes) {
	Con.df$cloneType <- NA
	for (x in seq_along(cloneTypes)) {
		names(cloneTypes)[x] <- paste0(
			names(cloneTypes[x]), ' (', cloneTypes[x-1], ' < X <= ',
			cloneTypes[x], ')'
		)
	}
	for (i in 2:length(cloneTypes)) {
		Con.df$cloneType <- ifelse(
			Con.df$Frequency>cloneTypes[i-1] & Con.df$Frequency<=cloneTypes[i],
			names(cloneTypes[i]),
			Con.df$cloneType
		)
	}
	Con.df
}

#' @title
#' Add scRepertoire-compatible clonotype information to a seurat object
#'
#' @description
#' This is a modified version of `scRepertoire::combineExpression`, which adds
#' the immune receptor information to a Seurat object's meta data. This function
#' also calculates the combined clonotype frequencies.
#'
#' Importantly, before using this function, ensure the barcodes of the seurat
#' object match the barcodes in the output of the
#' `scRepertoire::combinedContig()` call. Check `AddSampleBarcode()`
#' to change the prefixes of the Seurat object. If combining more than one
#' immune receptor type, barcodes with both receptors will be removed during the
#' combination process.
#'
#' @details
#' This function is a sort of "subset" of the original `scRepertoire` function
#' that combines and processes the clonotype data so that it can be used by
#' `APackOfTheClones` further downstream. The modified output seurat object
#' should be completely compatible to be used with other `scRepertoire`
#' functions. `scRepertoire::combineExpression` can also be ran again on the
#' modified seurat object and it shouldn't overwrite any of the existing data.
#'
#' @param contigs The product of `scRepertoire::CombineTCR()`,
#' `scRepertoire::CombineBCR()` or a list of both. (see scRepertoire's vignette)
#' @param seurat_obj The seurat object to attach the contigs to
#' @param cloneCall character. How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt) CDR3 amino acid (aa), or VDJC gene + CDR3 nucleotide
#' (strict).
#' @param chain character. Indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param cloneTypes named numeric. The bins for the grouping based on frequency
#' @param filterNA logical. Method to subset seurat object of barcodes without
#' clonotype information
#'
#' @importFrom dplyr %>%
#' @importFrom  rlang %||%
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
		contigs,
		seurat_obj,
		cloneCall = "strict",
		chain = "both",
		filterNA = FALSE,
		cloneTypes = c(
			Rare = 1e-4,Small = 0.001,Medium = 0.01,Large = 0.1,Hyperexpanded = 1
		)
) {
	call_time <- Sys.time()
	ce_warn_error_helper()

	options( dplyr.summarise.inform = FALSE )
	cloneTypes <- c(None = 0, cloneTypes)
	#df <- checkList(df)
	cloneCall <- clonecall_to_colname(cloneCall)
	Con.df <- NULL
	cell.names <- rownames(seurat_obj@meta.data)

	# first add the Frequency.all column for apotc
	data <- data.frame(dplyr::bind_rows(contigs), stringsAsFactors = FALSE)
	data <- data %>%
		dplyr::group_by(data[,cloneCall]) %>%
		dplyr::mutate(Frequency = dplyr::n())

	Con.df <- data[,c("barcode", "CTgene", "CTnt",
					  "CTaa", "CTstrict", "Frequency.all")]

	Con.df <- add_clonetype_porportion(Con.df, cloneTypes)
	PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt",
								"CTaa", "CTstrict", "Frequency", "cloneType")])
	dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
	PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
	rownames(PreMeta) <- PreMeta$barcode

	if (length(which(rownames(PreMeta) %in%
					 rownames(seurat_obj[[]])))/length(rownames(seurat_obj[[]])) < 0.01){
		warning(paste(
			"< 1% of barcodes match: Ensure the barcodes in the Seurat object",
			"match the barcodes in the combined immune receptor list from",
			"scRepertoire - most common issue is the addition of the prefixes",
			"corresponding to `samples` and 'ID' in the",
			"`scRepertoire::combineTCR`/`scRepertoire::BCR()` functions"
		))
	}
	seurat_obj[[names(PreMeta) %||% colnames(PreMeta)]] <- PreMeta

	if (filterNA) { seurat_obj <- filteringNA(seurat_obj, cloneTypes) }
	seurat_obj$cloneType <- factor( # dollar sign is for accessing metadata!
		seurat_obj$cloneType, levels = rev(names(cloneTypes))
	)
	# remove barcode col an the
	seurat_obj@meta.data[["barcode"]] <- NULL

	seurat_obj@commands[["combineSeuratExpression"]] <- make_apotc_cmd(
		call_time, seurat_obj@active.assay
	)

	seurat_obj
}

isCompatibleWithScrep <- function(combined_pbmc) {
	print("unfinished")
}
