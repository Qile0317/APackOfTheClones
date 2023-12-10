# script to manage the interface for accessing the apotc data
# could also put it in ApotcData.R

# from the input of RunAPOTC, convert the condition to a call to be put in
# @meta.data %>% dpylr::filter(eval(parse(text = "output of this func")))
# assume that metadata_filter is a valid ADDITIONAL filter condition.
#
# assume that varargs_list is a valid named list where each name is a column
# and element is a string vector of which factors to INCLUDE
parse_to_metadata_filter_str <- function(metadata_filter, varargs_list) {

    if (is.null(metadata_filter) && identical(varargs_list, list(NULL))) {
        return("")
    }

    filter_string <- ""
    colnames <- names(varargs_list)

    for (i in seq_along(varargs_list)) {
        
        sub_filter_cond_string <- col_cond_vec_to_filter_str(
            condition_vector = varargs_list[[i]], colname = colnames[i]
        )

        filter_string <- paste(
            filter_string, "(", sub_filter_cond_string, ") & ", sep = ""
        )
    }

    filter_string <- substr(filter_string, 1, length(filter_string) - 2)

    if (!is.null(metadata_filter)) {
        filter_string <- paste(
            "((", filter_string, ") & (", metadata_filter, "))", sep = ""
        )
    }

    filter_string
}

col_cond_vec_to_filter_str <- function(condition_vector, colname) {
  UseMethod("col_condition_vec_to_filter_string")
}

col_cond_vec_to_filter_str.character <- function(
    condition_vector, colname
) {
    col_condition_vec_to_filter_string_with_insert(
        condition_vector = condition_vector, colname = colname,
        insert_char = "'"
    )
}

col_cond_vec_to_filter_str.default <- function(
    condition_vector, colname
) {
    col_condition_vec_to_filter_string_with_insert(
        condition_vector = condition_vector, colname = colname,
        insert_char = ""
    )
}

col_condition_vec_to_filter_string_with_insert <- function(
    condition_vector, colname, insert_char
) {
    filter_str <- ""
    for (i in seq_along(condition_vector)) {
        filter_str <- paste(
            filter_str, "(",
            colname, " == ", insert_char, condition_vector[i], insert_char,
            ") | ", sep = ""
        )
    }
    substr(filter_str, 1, length(filter_str) - 2)
}

.defaultApotcDataObjId <- "__all__"
utils::globalVariables(c(".defaultApotcDataSample"))

# from the input of RunAPOTC, convert the condition to the apotc data sample id where
# its stored under @misc[["APackOfTheClones"]][[id]]
parse_to_object_id <- function(
    reduction_base, clonecall, varargs_list, metadata_filter
) {
	object_id <- paste("|", reduction_base, "|", clonecall, "|", sep = "")

    if (identical(varargs_list, list(NULL))) {
        object_id <- paste(object_id, "-|", sep = "")
    } else {
        object_id <- paste(
            object_id, varargs_list_to_id_segment(varargs_list), "|", sep = ""
        )
    }

    if (is.null(metadata_filter)) {
        return(paste(object_id, "-|", sep = ""))
    }
    paste(object_id, metadata_filter)
}

varargs_list_to_id_segment <- function(varargs_list) {
    # TODO
}

# make user getters for apotcData

# new format, there will be a list of apotc objects in the seurat@misc slot. the list will be named apotc.
# each is dependent on reduction/samples and within the list there will be named elements for each reduction/sample combo
# and make it optional in RunAPOTC if this should be stored. APOTCPlot will then be able to have the apotc obj slot input
# alternatively the sample/id configuration.
