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
#
# assume that varargs_list is a valid named list where each name is a column
# and element is a string vector of which factors to INCLUDE
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

parse_to_object_id <- function(
    reduction_base, clonecall, varargs_list, metadata_filter,
    .idSepStr = ";", .idNullStr = "_"
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

# getting and setting

containsApotcRun <- function(seurat_obj, obj_id) {
    return(!is.null(getApotcData(seurat_obj, obj_id)))
}

getApotcData <- function(seurat_obj, obj_id) {
    seurat_obj@misc[["APackOfTheClones"]][[obj_id]]
}

setApotcData <- function(seurat_obj, obj_id, apotc_obj) {
    seurat_obj@misc[["APackOfTheClones"]][[obj_id]] <- apotc_obj
    seurat_obj
}

#' @title
#' Get all object ids of previous RunAPOTC runs on a seurat object
#' 
#' @description
#' A convenience function to get all object ids of previous RunAPOTC run IDs
#' 
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#' 
#' @return a character vector of all object ids of previous RunAPOTC runs
#' @export
#' 
#' @example
#' # TODO
#' 
getApotcDataId <- function(seurat_obj) {
    if (!is_seurat_object(seurat_obj)) {
        stop("input must be a seurat object")
    }
    obj_list <- seurat_obj@misc[["APackOfTheClones"]]
    if (!isnt_empty(obj_list) && !is.null(obj_list)) {
        stop("No APackOfTheClones data found in seurat object")
    }
    names(obj_list)
}

#' @title
#' Get the object id of the last RunAPOTC run on a seurat object
#' 
#' @description
#' A convenience function to get the last object ids of the previous RunAPOTC
#' run, to be used by [APOTCPlot] and [AdjustAPOTC]
#' 
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#' 
#' @return a character of the object id of the last RunAPOTC call
#' @export
#' 
#' @example
#' # TODO
#' 
getLastApotcDataId <- function(seurat_obj) {
    getlast(getApotcDataId(seurat_obj))
}
