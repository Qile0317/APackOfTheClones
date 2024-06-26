#' @title
#' Check for the existence of an APackOfTheClones run with its run id
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A convenience function to check for the existence of an APackOfTheClones
#' run with its run id, regardless of if any run has been made
#'
#' @param seurat_obj a seurat object
#' @param run_id character. The id of the associated ApotcRun.
#'
#' @return A logical indicating whether the run exists.
#' @export
#'
#' @examples
#' pbmc <- RunAPOTC(
#'     seurat_obj = get(data("combined_pbmc")),
#'     reduction_base = "umap",
#'     clonecall = "strict",
#'     run_id = "run1",
#'     verbose = FALSE
#' )
#'
#' containsApotcRun(pbmc, "run1")
#' #> [1] TRUE
#'
#' containsApotcRun(pbmc, "run2")
#' #> [1] FALSE
#'
containsApotcRun <- function(seurat_obj, run_id) {

    if (!is_seurat_object(seurat_obj)) stop("input must be a seurat object")
    typecheck(run_id, is_a_character)

    if (!containsAnyApotcData(seurat_obj)) {
        return(FALSE)
    }

    any(getApotcDataIds(seurat_obj) == run_id)
}

#' @title
#' Rename an APackOfTheClones run
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A function to rename an APackOfTheClones run identified by its run id in
#' a Seurat object.
#'
#' @param seurat_obj A Seurat object containing APackOfTheClones data - the
#' output of [RunAPOTC]
#' @param old_run_id Character. The current id of the APackOfTheClones run to
#' be renamed.
#' @param new_run_id Character. The new id to assign to the APackOfTheClones
#' run.
#'
#' @return A Seurat object with the APackOfTheClones run renamed.
#' @export
#'
#' @examples
#' pbmc <- RunAPOTC(
#'     seurat_obj = get(data("combined_pbmc")),
#'     reduction_base = "umap",
#'     clonecall = "strict",
#'     run_id = "run1",
#'     verbose = FALSE
#' )
#'
#' pbmc <- renameApotcRun(pbmc, "run1", "new_run")
#' # Now "run1" has been renamed to "new_run"
#'
renameApotcRun <- function(seurat_obj, old_run_id, new_run_id) {

    if (!is_seurat_object(seurat_obj)) stop("input must be a seurat object")
    typecheck(old_run_id, is_a_character)
    typecheck(new_run_id, is_a_character)

    if (!containsApotcRun(seurat_obj, old_run_id))
        stop("There's no run named `", old_run_id, "`")
    if (containsApotcRun(seurat_obj, new_run_id))
        stop("There is already a run named `", new_run_id, "`")

    names(seurat_obj@misc$APackOfTheClones)[
        getApotcDataIds(seurat_obj) == old_run_id
    ] <- new_run_id

    seurat_obj
}

#' @title
#' Delete the results of an APackOfTheClones run
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A convenience function to erase all data associated with a particular run,
#' including the ApotcData and command in seurat_obj@command. The run_id would
#' be no longer accessible afterwards.
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
#'     seurat_obj = get(data("combined_pbmc")),
#'     reduction_base = "umap",
#'     clonecall = "strict",
#'     run_id = "run1",
#'     verbose = FALSE
#' )
#'
#' getApotcDataIds(pbmc)
#' #> [1] "run1"
#'
#' # delete the data
#' pbmc <- deleteApotcData(pbmc, "run1")
#'
#' getApotcDataIds(pbmc)
#' #> NULL
#'
deleteApotcData <- function(seurat_obj, run_id) {

    if (!is_seurat_object(seurat_obj)) stop("input must be a seurat object")
    if (length(run_id) != 1) stop("the `run_id` argument must be of length 1")
    if (!containsApotcRun(seurat_obj, run_id)) {
        stop(paste("no run with id:", run_id, "is present"))
    }

    seurat_obj <- setApotcData(seurat_obj, run_id, NULL)
    seurat_obj@commands[[get_command_name("RunAPOTC", run_id)]] <- NULL
    seurat_obj
}

#' @title
#' Get all run ids of previous RunAPOTC runs on a seurat object
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A convenience function to get all run ids of previous RunAPOTC run IDs
#'
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#'
#' @return a character vector of all run ids of previous RunAPOTC runs, in
#' the order they were ran in. If there are no runs on the object, it returns
#' `NULL`.
#' @export
#'
#' @examples
#' pbmc <- RunAPOTC(
#'     seurat_obj = get(data("combined_pbmc")),
#'     reduction_base = "umap",
#'     clonecall = "strict",
#'     verbose = FALSE
#' )
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
    ids <- names(seurat_obj@misc[["APackOfTheClones"]])
    if (identical(ids, character(0)) || is.null(ids)) return(NULL)
    ids
}

#' @title
#' Get the object id of the most recent RunAPOTC run on a seurat object
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A convenience function to get the object id of the most recent valid
#' [RunAPOTC] run, to be used by [APOTCPlot] and [AdjustAPOTC]
#'
#' @param seurat_obj a seurat object that has had RunAPOTC ran on it before in
#' order of the functions being called.
#'
#' @return a character of the object id of the last [RunAPOTC] call
#' @export
#'
#' @examples
#' # first run
#' pbmc <- RunAPOTC(
#'     seurat_obj = get(data("combined_pbmc")),
#'     reduction_base = "umap",
#'     clonecall = "strict",
#'     verbose = FALSE
#' )
#'
#' getLastApotcDataId(pbmc)
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
#' getLastApotcDataId(pbmc)
#' #> [1] "umap;CTgene;_;_"
#'
getLastApotcDataId <- function(seurat_obj) {
    getlast(getApotcDataIds(seurat_obj))
}

# TODO change object ids-make it all nonempty arguments

# object id generation:

# from the input of RunAPOTC, convert the condition to a call to be put in
# @meta.data %>% dpylr::filter(eval(parse(text = "output of this func")))
# assume that metadata_filter is a valid ADDITIONAL filter condition.
# assume that varargs_list is a valid named list where each name is a column
# and element is a string vector of which factors to INCLUDE
# In hindsight %in% can be used but its now unessecary to make changes
parse_to_metadata_filter_str <- function(metadata_filter, varargs_list) {

    if (!is_valid_args(varargs_list)) {
        if (is_valid_filter_str(metadata_filter)) {
            return(strip_unquoted_spaces(metadata_filter))
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

    # TODO FIXME - this will allow symbolic filtering as well
    #metadata_filter <- convert_to_string_if_expr(deparse(substitute(metadata_filter)))

    if (is_valid_filter_str(metadata_filter)) {
        filter_string <- paste(
            "(", filter_string, ")&(", metadata_filter, ")", sep = ""
        )
    }

    strip_unquoted_spaces(filter_string)
}

is_valid_args  <- function(varargs_list) isnt_empty(varargs_list)

is_valid_filter_str <- function(metadata_string) {
    if (is.null(metadata_string)) return(FALSE)
    if (identical(strip_spaces(metadata_string), "")) return(FALSE)
    return(TRUE)
}

col_cond_vec_to_filter_str <- function(condition_vector, colname) {

    insertchar <- ""
    if (is.character(condition_vector)) insertchar <- "'"

    col_conds_to_str_w_insert(
        condition_vector = condition_vector,
        colname = colname,
        insert_char = insertchar
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

# if symbolic, convert to filter string
convert_to_string_if_expr <- function(deparsed_metadata_filter) {
   gsub("\"", "'", deparsed_metadata_filter)
}

# functions for converting args of RunAPOTC to the apotc data sample id
# stored under under @misc[["APackOfTheClones"]][[id]]
infer_object_id_if_needed <- function(args, varargs_list) {
    if (!should_compute(args$run_id)) return(args$run_id)
    
    if (
        is.null(args$reduction_base) &&
            is.null(args$clonecall) &&
            is.null(args$extra_filter) &&
            is_empty(varargs_list)
    ) {
        latest_id <- getLastApotcDataId(args$seurat_obj)

        # verbosity?
        message(paste(
            "* using the latest APackOfTheClones Run Data with object id:",
            latest_id
        ))
        return(latest_id)
    }

    parse_to_object_id(
        reduction_base = attempt_correction(
            args$seurat_obj, args$reduction_base
        ),
        clonecall = .theCall(args$seurat_obj@meta.data, args$clonecall),
        varargs_list = varargs_list,
        metadata_filter = args$extra_filter
    )
}

utils::globalVariables(c(".idSepStr", ".idNullStr"))
.idSepStr = ";" # not the best idea - `;` on its own is end of line.
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

# deparse into a list of
# reduction_base, clonecall, varargs_list, metadata_filter
deparse_object_id <- function(object_id) {

    segments <- strsplit(object_id, .idSepStr)[[1]]

    reduction_base <- segments[1]
    clonecall <- segments[2]

    if (identical(segments[3], .idNullStr)) {
        varargs_list <- list()
    } else {
        varargs_list <- id_segment_to_varargs_list(segments[3])
    }

    if (identical(segments[4], .idNullStr)) {
        metadata_filter <- ""
    } else {
        metadata_filter <- segments[4]
    }

    list(
        reduction_base = reduction_base, clonecall = clonecall,
        varargs_list = varargs_list, extra_filter = metadata_filter
    )
}

id_segment_to_varargs_list <- function(id_segment) {
    eval(as_expression("list(", id_segment, ")"))
}

# getting and setting related functions

getApotcData <- function(seurat_obj, obj_id) {
    seurat_obj@misc[["APackOfTheClones"]][[obj_id]]
}

getLastApotcData <- function(seurat_obj) {
    getApotcData(seurat_obj, getLastApotcDataId(seurat_obj))
}

getApotcDataIfExistsElseCreate <- function(seurat_obj, obj_id, arglist, ...) {

    if (!is.null(obj_id)) {
        if (!containsApotcRun(seurat_obj, obj_id)) {
            stop(call. = FALSE, "No runs found with the id: `", obj_id, "`")
        }
        return(getApotcData(seurat_obj, obj_id))
    }

    if (is_an_integer(arglist)) arglist <- parent.frame(arglist)
    if (is.null(arglist)) arglist <- parent.frame()

    actual_metadata_filter <- parse_to_metadata_filter_str(
        arglist$extra_filter, list(...)
    )

    seurat_obj %>%
        RunAPOTC(
            reduction_base = arglist$reduction_base,
            clonecall = arglist$clonecall,
            extra_filter = actual_metadata_filter,
            alt_ident = arglist$alt_ident,
            run_id = "__getApotcDataIfExistsElseCreate__",
            override = TRUE,
            verbose = FALSE
        ) %>%
        getApotcData("__getApotcDataIfExistsElseCreate__")
}

setApotcData <- function(seurat_obj, obj_id, apotc_obj) {
    seurat_obj@misc[["APackOfTheClones"]][[obj_id]] <- apotc_obj
    seurat_obj
}

containsAnyApotcData <- function(seurat_obj) {
    !is.null(getApotcDataIds(seurat_obj))
}
