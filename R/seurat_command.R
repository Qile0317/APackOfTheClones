# script for the function to add the APOTC command records to seurat objects
# Is copied from SeuratObject R/command.R script, under the MIT license

# almost identical to seurat version except only data.frame's names are saved in
# params and assay.used is only "RNA"
make_apotc_command <- function(call_time = "auto", assay = "RNA") {
    if (should_assume(call_time)) {call_time <- Sys.time()}

    if (as.character(x = sys.calls()[[1]])[1] == "do.call") {
        call_string <- deparse(expr = sys.calls()[[1]])
        command_name <- as.character(x = sys.calls()[[1]])[2]
    } else {
        command_name <- as.character(
            x = deparse(expr = sys.calls()[[sys.nframe() - 1]])
        )
        command_name <- gsub(
            pattern = "\\.Seurat",
            replacement = "",
            x = command_name
        )
        call_string <- command_name
        command_name <- seurat_extractfield(
            string = command_name,
            field = 1,
            delim = "\\("
        )
    }

    argnames <- names(formals(fun = sys.function(which = sys.parent(n = 1))))
    argnames <- process_argnames(argnames)
    params <- get_parent_params(n = 2)

    # return the command object
    methods::new(
        Class = 'SeuratCommand',
        name = command_name,
        params = params,
        time.stamp = call_time,
        call.string = call_string,
        assay.used = assay
    )
}

seurat_extractfield <- function(string, field = 1, delim = "_") {
    fields <- as.numeric(
        x = unlist(x = strsplit(x = as.character(x = field), split = ","))
    )
    if (length(x = fields) == 1) {
        return(strsplit(x = string, split = delim)[[1]][field])
    }
    return(paste(
        strsplit(x = string, split = delim)[[1]][fields], collapse = delim
    ))
}

# function to be used within another parent function, extracting the arguments
# to the parent function and returning it as a named list, while also allowing
# filtering out of certain object types to save memory

process_argnames <- function(argnames) {
    argnames <- BiocGenerics::grep(
        pattern = "object",
        x = argnames,
        invert = TRUE,
        value = TRUE
    )
    argnames <- BiocGenerics::grep(
        pattern = "anchorset",
        x = argnames,
        invert = TRUE,
        value = TRUE
    )
    argnames <- BiocGenerics::grep(
        pattern = "\\.\\.\\.",
        x = argnames,
        invert = TRUE,
        value = TRUE
    )
    argnames
}

get_parent_params <- function(
    n = 1,
    excluded_types = c("Seurat", "list"),
    only_named_types = c("data.frame", "data.table")
) {
    argnames <- names(formals(fun = sys.function(which = sys.parent(n = n))))
    argnames <- process_argnames(argnames)

    params <- list()
    p.env <- parent.frame(n = n)
    argnames <- intersect(x = argnames, y = ls(name = p.env))
    for (arg in argnames) {
        param_value <- get(x = arg, envir = p.env)

        is_excluded_type <- FALSE
        for (obj_type in excluded_types) {
            if (inherits(param_value, obj_type)) {
                is_excluded_type <- TRUE
                break
            }
        }
        if (is_excluded_type) {next}

        for (obj_type in only_named_types) {
            if (inherits(param_value, obj_type)) {
                param_value <- names(param_value)
                break
            }
        }
        params[[arg]] <- param_value
    }
    params
}

# convenience function (should remove)
get_cmd <- function(seurat_obj, cmd) {
    seurat_obj@commands[["RunAPOTC"]]@params[[cmd]]
}
