# script for the function to add the APOTC command records to seurat objects
# Is copied from SeuratObject R/command.R script, under the MIT license

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

# almost identical to seurat version except only data.frame's names are saved in
# params and assay.used is only "RNA"
make_apotc_command <- function(call_time = "auto") {

    if (identical(call_time, "auto")) {call_time <- Sys.time()}

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

    argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
    argnames <- grep(
        pattern = "object",
        x = argnames,
        invert = TRUE,
        value = TRUE
    )
    argnames <- grep(
        pattern = "anchorset",
        x = argnames,
        invert = TRUE,
        value = TRUE
    )
    argnames <- grep(
        pattern = "\\.\\.\\.",
        x = argnames,
        invert = TRUE,
        value = TRUE
    )

    params <- list()
    p.env <- parent.frame(n = 1)
    argnames <- intersect(x = argnames, y = ls(name = p.env))
    for (arg in argnames) {
        param_value <- get(x = arg, envir = p.env)
        if (inherits(x = param_value, what = 'Seurat')) {
            next
        }
        if (inherits(x = param_value, what = 'data.frame')) {
            params[[arg]] <- names(param_value)
            next
        }
        params[[arg]] <- param_value
    }

    command_name <- sub(
        pattern = "[\\.]+$",
        replacement = "",
        x = command_name,
        perl = TRUE
    )
    command_name <- sub(
        pattern = "\\.\\.", replacement = "\\.", x = command_name, perl = TRUE
    )

    # return the command object
    methods::new(
        Class = 'SeuratCommand',
        name = command_name,
        params = params,
        time.stamp = call_time,
        call.string = call_string,
        assay.used = "RNA"
    )
}

# convinience function
get_cmd <- function(seurat_obj, cmd) {
    seurat_obj@commands[["RunAPOTC"]]@params[[cmd]]
}
