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

# convinience function (should remove)
get_cmd <- function(seurat_obj, cmd) {
    seurat_obj@commands[["RunAPOTC"]]@params[[cmd]]
}
