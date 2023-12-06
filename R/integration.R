ce_warn_error_helper <- function() {
    params <- get_parent_params(n = 2)
    print(params) # unfinished, this is just a test
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

    if (x %in% c("CTnt", "CTgene", "CTaa", "CTstrict")) {
        return(x)
    }
    x <- tolower(x)
    if (!is.null(.clonecall_dictionary[[x]])) {
        return(.clonecall_dictionary[[x]])
    }
    stop(paste(
        "input cloneCall invalid, did you mean: '",
        closest_word(x, names(.clonecall_dictionary)),
        "'?",
        sep = ""
    ))
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
#' @param filterNA logical. Method to subset seurat object of barcodes without
#' clonotype information
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
    filterNA = FALSE
) {
    call_time <- Sys.time()
    ce_warn_error_helper()
    cloneCall <- clonecall_to_colname(cloneCall)

    is_first_call <- TRUE
    if (!is.null(seurat_obj@commands[["combineExpression"]])) {
        is_first_call <- FALSE
    }

    seurat_obj <- scRepertoire_combineExpression(
        df = contigs, sc = seurat_obj, cloneCall = cloneCall, chain = chain,
        group.by = NULL, filterNA = filterNA, addLabel = TRUE
    )

    seurat_obj@meta.data[["barcode"]] <- NULL
    if (is_first_call) {
        seurat_obj@commands[["combineExpression"]] <- NULL
    }
    seurat_obj@commands[["combineSeuratExpression"]] <- make_apotc_cmd(
        call_time, seurat_obj@active.assay
    )

    seurat_obj
}
