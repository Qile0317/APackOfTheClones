AddSamplePrefix_warner <- function(sc_list, params, sep) {
	if (identical(sep, "")) {
		stop("the separator `sep` probably should not be ''")
	}
	if (!isnt_empty(sc_list)) {
		stop("sc_list has to be a list of valid objects!")
	}
	num_samples <- length(sc_list)
	for (i in seq_along(params)) {
		if (!is.character(params[[i]])) {
			stop("every parameter has to be a character vector")
		}
		if (length(params[[i]]) != num_samples) {
			stop(paste(
				"every parameter has to be the same length as the number of",
				"seurat objects"
			))
		}
	}
}

get_barcode_prefix_end_index <- function(sc_list, sep = "_") {
	last_occurence_index(get_rna_assay_barcodes(sc_list[[1]])[1], sep)
}

construct_param_list <- function(samples, ID, extra_params) {
	params <- list()
	if (!is.null(samples)) {
		params <- c(params, list("sample" = samples))
	}
	if (!is.null(ID)) {
		params <- c(params, list("ID" = ID))
	}
	if (!is.null(extra_params)) {
		params <- c(params, extra_params)
	}
	params
}

add_repped_2d_list_row_to_df <- function(
	df, l, row_index, err_msg = paste("Some or all of the parameter names",
									  "already exist in one or all of the",
									  "seurat object(s)")
) {
	list_names <- names(l)
	objs_to_be_added <- extract_2d_list_row(l, row_index)
	df_height <- nrow(df)

	if (has_common_strs(colnames(df), list_names)) {stop(err_msg)}
	for (i in seq_along(l)) {
		df[[list_names[i]]] <- rep(objs_to_be_added[i], df_height)
	}
	df
}

add_params_to_metadata <- function(seurat_obj, params, ind) {
	seurat_obj@meta.data <- add_repped_2d_list_row_to_df(
		seurat_obj@meta.data, params, ind
	)
}

change_metadata_barcode <- function(seurat_obj, new_barcode) {
	rownames(seurat_obj@meta.data) <- new_barcode
	seurat_obj
}

change_assay_barcode <- function(seurat_obj, new_barcode) {
	for (i in 1:seq_along(seurat_obj@assays)) {
		if (!inherits(seurat_obj@assays[[i]], "Assay")) {
			next
		}
		try(seurat_obj@assays[[i]]@data@Dimnames[[2]] <- new_barcode, TRUE)
	}
	seurat_obj
}

change_reduction_barcode <- function(seurat_obj, new_barcode) {
	for (i in 1:seq_along(seurat_obj@reductions)) {
		if (!inherits(seurat_obj@reductions[[i]], "DimReduc")) {
			next
		}
		try(rownames(seurat_obj@reductions[[i]]) <- new_barcode, TRUE)
	}
	seurat_obj
}

change_barcode <- function(seurat_obj, new_barcode) {
	seurat_obj <- change_metadata_barcode(seurat_obj, new_barcode)
	seurat_obj <- change_assay_barcode(seurat_obj, new_barcode)
	seurat_obj <- change_reduction_barcode(seurat_obj, new_barcode)
	seurat_obj
}

#' @title
#' Add scRepertoire prefixes to barcodes of a seurat object(s) pre-integration
#'
#' @description
#' placeholder.
#'
#' @param sc_list A `list` of seurat objects to be modified
#' @param samples unfinished
#' @param ID unfinished
#' @param ... unfinished
#' @param add_to_metadata logical. If `TRUE`, adds the entire processed clonotype data
#' by barcode into the seurat object meta.data dataframe.
#' @param sep character. Decides the string to join new sample level informations.
#' Defaults to `"_"`.
#'
#' @return
#' The modified version of the input list of seurat objects `sc_list`
#'
#' @examples
#' data("mini_seurat_obj", "")
#'
#' @export
#'
AddSamplePrefix <- function(
	sc_list, samples = NULL, ID = NULL, ... = NULL, add_to_metadata = FALSE,
	sep = "_"
) {
	params <- construct_param_list(samples, ID, list(...))
	AddSamplePrefix_warner(sc_list, params, sep)
	prefix_vector <- construct_prefix_vector(params, sep)

	# adding prefixes and renaming list
	names(sc_list) <- prefix_vector
	for (i in seq_along(sc_list)) {
		new_barcodes <- paste0(
			prefix_vector[i], sep, rownames(sc_list[[i]]@meta.data)
		)

		sc_list[[i]] <- change_barcode(sc_list[[i]], new_barcodes)

		if (!add_to_metadata) {
			sc_list[[i]] <- add_params_to_metadata(sc_list[[i]], params, i)
		}
	}
	sc_list
}

#' @title
#' Remove barcode prefixes from a seurat object list
#'
#' @description
#' placeholder
#'
#' @param sc_list A `list` of seurat objects
#' @param sep placeholder
#'
#' @return
#' The modified version of the input list of seurat objects `sc_list`
#'
#' @examples
#' data("mini_seurat_obj", "")
#'
#' @export
#'
RemoveSamplePrefix <- function(sc_list, prefix_end = "auto", sep = "_") {
	if (should_assume(prefix_end)) {
		barcode_prefix_end_index <- get_barcode_prefix_end_index(sc_list, sep)
	} else if (is.numeric(prefix_end)) {
		barcode_prefix_end_index <- prefix_end
	}

	for (i in seq_along(sc_list)) {
		new_barcodes <- substring(
			names(sc_list[[i]]@meta.data), barcode_prefix_end_index + 2
		)
		sc_list[[i]] <- change_barcode(sc_list[[i]], new_barcodes)
	}
	sc_list
}

#' @title
#' Modify existing barcode prefixes from a seurat object list
#'
#' @description
#' placeholder
#'
#' @param sc_list A `list` of seurat objects
#' @param sep placeholder
#'
#' @return
#' The modified version of the input list of seurat objects `sc_list`
#'
#' @examples
#' data("mini_seurat_obj", "")
#'
#' @export
#'
ChangeSamplePrefix <- function(
	sc_list, samples = NULL, ID = NULL, ... = NULL,  add_to_metadata = FALSE,
	original_sep = "_", sep = "_"
) {
	AddSamplePrefix(
		RemoveSamplePrefix(sc_list, original_sep), samples, ID, ...,
		add_to_metadata, sep
	)
}

# to not have old code in vain, will make utility function to add in the raw
# contig list here. basically vectorized integrate_tcr but also need it to work
# for integrate bcr. could also make S3 dispatched methods for seurat and sc obj
