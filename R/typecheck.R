# TODO could make a currying function for length checking

# this function assumes its being called in an error checker
# and the error checker has an argument varargs_list
get_parent_func_args <- function(dn = 0L) {
    parent.frame(2L + dn) %>%
        as.list() %>%
        append(parent.frame()$varargs_list) %>%
        hash::hash()
}

is_pair <- function(x, type_checker) {
    if (length(x) != 2) return(FALSE)
    all(sapply(x, type_checker))
}

is_vector <- function(x) is.character(x) || is.numeric(x) || is.logical(x)

check_is_list_and_elements <- function(
    x, elem_type_checker, num_elements = NULL
) {
    if (!is.list(x)) return(FALSE)
    if (!is.null(num_elements)) if (length(x) != num_elements) return(FALSE)
    all(sapply(x, elem_type_checker))
}

is_seurat_object <- function(obj) inherits(obj, "Seurat")

is_an_apotc_ggplot <- isApotcGGPlot

is_a_character <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.character(x)
}

is_character <- is.character

is_a_logical <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.logical(x)
}

is_a_numeric <- function(x) {
    if (length(x) != 1) return(FALSE)
    is.numeric(x)
}

is_numeric_pair <- function(x) is_pair(x, is_a_numeric) && is.numeric(x)

is_list_of_numeric_pair <- function(x) {
    check_is_list_and_elements(x, is_numeric_pair)
}

is_a_positive_numeric <- function(x) {
    if (!is_a_numeric(x)) return(FALSE)
    x > 0
}

is_positive_numeric <- function(x) {
    if (!is_vector(x)) return(FALSE)
    if (length(x) < 1L) return(FALSE)
    all(sapply(x, is_a_positive_numeric))
}

is_an_integer <- function(x) {
    if (identical(x, Inf) || identical(x, -Inf)) return(TRUE)
    if (!is_a_numeric(x)) return(FALSE)
    as.numeric(x) == as.numeric(as.integer(x))
}

is_integer_pair <- function(x) is_pair(x, is_an_integer) && is.numeric(x)

is_integer <- function(x) {
    if (!is_vector(x)) return(FALSE)
    all(sapply(x, is_an_integer))
}

is_a_positive_integer <- function(x) {
    if (!is_an_integer(x)) return(FALSE)
    x > 0L
}

is_positive_integer <- function(x) {
    if (!is_integer(x)) return(FALSE)
    all(sapply(x, function(x) x > 0L))
}

# primary typechecking function where ... is a list of functions
typecheck <- function(x, ...) {

    typechecker_list <- list(...)
    if (any(sapply(typechecker_list, function(f) f(x)))) return()

    typechecker_str_vec <- convert_list_var_str_to_char(
        deparse(substitute(list(...)))
    )

    stop(call. = FALSE,
        "`",
        x %>% substitute() %>% deparse() %>% get_right_of_dollarsign(),
        "` ",
        create_err_msg(typechecker_str_vec)
    )
}

convert_list_var_str_to_char <- function(list_var_str) {
    substr(list_var_str, 6, nchar(list_var_str) - 1) %>%
        strsplit(", ") %>%
        getlast() %>%
        as.character()
}

create_err_msg <- function(typechecker_str_vec) {
    paste(
        "must be",
        get_error_strings(typechecker_str_vec) %>%
            join_error_strings()
    )
}

get_error_strings <- function(funcstrs) {
    sapply(
        funcstrs,
        function(s) prepend_indefinite_article(get_err_type_str(s))
    ) %>% sort()
}

get_err_type_str <- function(function_name_str) {

    if (function_name_str == "is.null") return("NULL")

    type <- gsub(
        "_", " ",
        sub("^is_(a|an)?_?", "", function_name_str)
    )

    if (grepl("^is_list_of_.*$", function_name_str)) {
        return(paste(type, "s", sep = ""))
    }

    if (grepl("^is_an?_.*", function_name_str)) {
        return(paste(type, "of length 1"))
    }

    if (grepl("^is_.*_pair$", function_name_str)) {
        return(type)
    }

    if (grepl("^is_.*$", function_name_str)) {
        return(paste(type, "vector"))
    }

    stop("dev error: pattern matching for typecheck failed")
}

join_error_strings <- function(error_string_vec) {

    num_error_strings <- length(error_string_vec)

    if (num_error_strings == 1) return(error_string_vec)

    if (num_error_strings == 2) {
        return(paste(error_string_vec, collapse = ", or "))
    }

    paste(
        paste(error_string_vec[1:num_error_strings - 1], collapse = ", "),
        getlast(error_string_vec), sep = ", or "
    )
}
