# script to manage the interface for accessing the apotc data
# all functions assume arguments are correct

is_valid_filter_str <- function(metadata_string) {
    if (is.null(metadata_string)) return(FALSE)
    if (identical(strip_spaces(metadata_string), "")) return(FALSE)
    return(TRUE)
}

is_valid_args  <- function(varargs_list) {
    isnt_empty(varargs_list)
}

# from the input of RunAPOTC, convert the condition to a call to be put in
# @meta.data %>% dpylr::filter(eval(parse(text = "output of this func")))
# assume that metadata_filter is a valid ADDITIONAL filter condition.
# assume that varargs_list is a valid named list where each name is a column
# and element is a string vector of which factors to INCLUDE
#
parse_to_metadata_filter_str <- function(metadata_filter, varargs_list) {

    if (!is_valid_args(varargs_list)) {
        if (is_valid_filter_str(metadata_filter)) {
            return(strip_spaces(metadata_filter))
        }
        return("")
    }

    filter_strings <- vector("character", length(varargs_list))
    colnames <- names(varargs_list)

    for (i in seq_along(varargs_list)) {
        filter_strings[i] <- col_cond_vec_to_filter_str(
            condition_vector = sort(unique(varargs_list[[i]])), colnames[i]
        )
    }

    filter_string <- sort_and_join_conds_by_and(filter_strings)

    if (is_valid_filter_str(metadata_filter)) {
        filter_string <- paste(
            "(", filter_string, ")&(", metadata_filter, ")", sep = ""
        )
    }

    strip_spaces(filter_string)
}

col_cond_vec_to_filter_str <- function(condition_vector, colname) {
  UseMethod("col_cond_vec_to_filter_str")
}

col_cond_vec_to_filter_str.character <- function(
    condition_vector, colname
) {
    col_conds_to_str_w_insert(
        condition_vector = condition_vector, colname = colname,
        insert_char = "'"
    )
}

col_cond_vec_to_filter_str.default <- function(
    condition_vector, colname
) {
    col_conds_to_str_w_insert(
        condition_vector = condition_vector, colname = colname,
        insert_char = ""
    )
}

col_conds_to_str_w_insert <- function(
    condition_vector, colname, insert_char
) {
    filter_str <- ""
    for (i in seq_along(condition_vector)) {
        filter_str <- paste(
            filter_str,
            colname, "==", insert_char, condition_vector[i], insert_char,
            "|", sep = ""
        )
    }
    substr(filter_str, 1, nchar(filter_str) - 1)
}

sort_and_join_conds_by_and <- function(filter_strings) {
    if (length(filter_strings) == 1) return(filter_strings[1])
    paste(
        "(", paste(sort(filter_strings), collapse = ")&("), ")", sep = ""
    )
}

# functions for converting args of RunAPOTC to the apotc data sample id
# stored under under @misc[["APackOfTheClones"]][[id]]

infer_object_id <- function(args, varargs_list) {
    if (
        is.null(args$reduction_base) &&
            is.null(args$clonecall) &&
            is.null(args$extra_filter) &&
            is_empty(varargs_list)
    ) {
        latest_id <- getLastApotcDataId(args$seurat_obj)
        message(paste(
            "* using the latest APackOfTheClones Run Data with object id:",
            latest_id
        ))
        return(latest_id)
    }

    parse_to_object_id(
        reduction_base = attempt_correction(args$reduction_base),
        clonecall = .theCall(args$seurat_obj@meta.data, args$clonecall),
        varargs_list = varargs_list,
        metadata_filter = args$extra_filter
    )
}

utils::globalVariables(c(".idSepStr", ".idNullStr"))
.idSepStr = ";"
.idNullStr = "_"

parse_to_object_id <- function(
    reduction_base, clonecall, varargs_list, metadata_filter
) {
	object_id <- paste(reduction_base, .idSepStr, clonecall, .idSepStr, sep = "")

    if (!is_valid_args(varargs_list)) {
        object_id <- paste(object_id, .idNullStr, .idSepStr, sep = "")
    } else {
        object_id <- paste(
            object_id, varargs_list_to_id_segment(varargs_list),
            .idSepStr, sep = ""
        )
    }

    if (!is_valid_filter_str(metadata_filter)) {
        return(paste(object_id, .idNullStr, sep = ""))
    }
    paste(object_id, gsub(" ", "", metadata_filter), sep = "")
}

get_default_apotc_id <- function(reduction_base, clonecall) {
    parse_to_object_id(
        reduction_base = reduction_base, clonecall = clonecall,
        varargs_list = list(), metadata_filter = NULL
    )
}

varargs_list_to_id_segment <- function(varargs_list) {
    segments <- vector("character", length(varargs_list))
    colnames <- names(varargs_list)
    for (i in seq_along(varargs_list)) {
        segments[i] <- paste(
            colnames[i], "=",
            repr_as_string(sort(unique(varargs_list[[i]]))),
            sep = ""
        )
    }

    if (length(segments) == 1) {
        return(segments)
    }
    paste(sort(segments), collapse = ",")
}

#obj_id_to_readable_str

# getting and setting related functions

containsApotcRun <- function(seurat_obj, obj_id) {
    return(!is.null(getApotcData(seurat_obj, obj_id)))
}

user_contains_apotc_run <- function(seurat_obj, obj_id) {
    # TODO
    containsApotcRun(obj_id)
}

getApotcData <- function(seurat_obj, obj_id) {
    seurat_obj@misc[["APackOfTheClones"]][[obj_id]]
}

getLastApotcData <- function(seurat_obj) {
    getApotcData(seurat_obj, getLastApotcDataId(seurat_obj))
}

setApotcData <- function(seurat_obj, obj_id, apotc_obj) {
    seurat_obj@misc[["APackOfTheClones"]][[obj_id]] <- apotc_obj
    seurat_obj
}

#' @title
#' Delete the results of an APackOfTheClones run
#' 
#' @description
#' A convenience function to erase all data associated with a particular run,
#' including the ApotcData and command in seurat_obj@command.
#' 
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#' @param run_id character. The id of the associated ApotcRun.
#' 
#' @return The modified input seurat object
#' @export
#' 
#' @examples
#' pbmc <- RunAPOTC(
# '     seurat_obj = get(data("combined_pbmc")),
# '     reduction_base = "umap",
# '     clonecall = "strict",
#'      run_id = "run1"
# '     verbose = FALSE
# ' )
#' 
#' getApotcDataIds(pbmc)
#' #> [1] "run1"
#' 
#' pbmc <- deleteApotcData(seurat_obj, run_id)
#' 
#' getApotcDataIds(pbmc)
#' #> NULL
#' 
deleteApotcData <- function(seurat_obj, run_id) {

    if (!is_seurat_object(seurat_obj)) stop("input must be a seurat object")
    if (length(run_id) != 1) stop("the `run_id` argument must be of length 1")
    user_contains_apotc_run(seurat_obj, run_id)

    seurat_obj <- setApotcData(seurat_obj, run_id, NULL)
    seurat_obj@commands[[get_command_name("RunAPOTC", run_id)]] <- NULL
    seurat_obj
}

#' @title
#' Get all run ids of previous RunAPOTC runs on a seurat object
#'
#' @description
#' A convenience function to get all run ids of previous RunAPOTC run IDs
#'
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#'
#' @return a character vector of all run ids of previous RunAPOTC runs, in
#' the order they were ran in.
#' @export
#'
#' @examples
# ' pbmc <- RunAPOTC(
# '     seurat_obj = get(data("combined_pbmc")),
# '     reduction_base = "umap",
# '     clonecall = "strict",
# '     verbose = FALSE
# ' )
#'
#' getApotcDataIds(pbmc)
#' #> [1] "umap;CTstrict;_;_"
#'
#' pbmc <- RunAPOTC(
#'     seurat_obj = pbmc,
#'     reduction_base = "umap",
#'     clonecall = "gene",
#'     verbose = FALSE
#' )
#'
#' getApotcDataIds(pbmc)
#' #> [1] "umap;CTstrict;_;_" "umap;CTgene;_;_"
#'
getApotcDataIds <- function(seurat_obj) {
    if (!is_seurat_object(seurat_obj)) stop("input must be a seurat object")
    obj_list <- seurat_obj@misc[["APackOfTheClones"]]
    if (!isnt_empty(obj_list) && !is.null(obj_list)) {
        stop("No APackOfTheClones data found in seurat object")
    }
    names(obj_list)
}

#' @title
#' Get the object id of the most recent RunAPOTC run on a seurat object
#'
#' @description
#' A convenience function to get the object id of the most recent valid
#' [RunAPOTC] run, to be used by [APOTCPlot] and [AdjustAPOTC]
#'
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#'
#' @return a character of the object id of the last [RunAPOTC] call
#' @export
#'
#' @example
#' # first run
# ' pbmc <- RunAPOTC(
# '     seurat_obj = get(data("combined_pbmc")),
# '     reduction_base = "umap",
# '     clonecall = "strict",
# '     verbose = FALSE
# ' )
#'
#' getApotcDataIds(pbmc)
#' #> [1] "umap;CTstrict;_;_"
#'
#' # second run with a different clonecall
#' pbmc <- RunAPOTC(
#'     seurat_obj = pbmc,
#'     reduction_base = "umap",
#'     clonecall = "gene",
#'     verbose = FALSE
#' )
#'
#' getApotcDataIds(pbmc)
#' #> [1] "umap;CTgene;_;_"
#'
getLastApotcDataId <- function(seurat_obj) {
    getlast(getApotcDataIds(seurat_obj))
}
