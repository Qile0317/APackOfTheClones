#' @title
#' Run the APackOfTheClones method on a combined Seurat object for
#' downstream visualization of clonal expansion
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Computes necessary information for an APackOfTheClones
#' clonal expansion plot ([APOTCPlot]) and stores it in the seurat object.
#' Gets sizes of unique clones and utilizes a circle-packing algorithm to
#' pack circles representing individual clones in approximately the same
#' dimensional reduction (`reduction_base`) coordinates.
#'
#' The parameter `extra_filter` along with an unlimited number of additional
#' keyword arguments can be used to filter the cells by certain conditions
#' in the metadata, and new results will be stored in addition to other runs
#' the users may have done.
#'
#' Each APackOfTheClones run is uniquely identified by the parameters
#' `reduction_base`, `clonecall`, `extra_filter`, and any additional keywords
#' passed to filter the metadata. Each distinct run result is stored in the
#' seurat object and has an associated Id generated from the aforementioned
#' parameters. To view the id of the latest run, call [getLastApotcDataId].
#' To view all the ids of previous runs, call [getApotcDataIds]. To work further
#' with a specific run (most importantly, plotting), the user can use this id
#' in the arguments with is slightly more convenient than passing in the
#' original RunAPOTC parameters again but both ways work.
#'
#' If the user wishes to manually customize/fix the expansion plot
#' generated, the circular packing information can be modified
#' with the [AdjustAPOTC] function.
#'
#' @param seurat_obj Seurat object with one or more dimension reductions and
#' already have been integrated with a TCR/BCR library with
#' [combineSeuratExpression] or [scRepertoire::combineExpression]
#' @param reduction_base character. The seurat reduction to base the clonal
#' expansion plotting on. Defaults to `'umap'` but can also be `'tsne'` or
#' `'pca'`. If `'pca'``, the cluster coordinates will be based on PC1 and PC2.
#' However, generally APackOfTheClones is used for displaying UMAP and
#' occasionally t-SNE versions to intuitively highlight clonal expansion.
#' @param clonecall character. The column name in the seurat object metadata to
#' use. See `scRepertoire` documentation for more information about this
#' parameter that is central to both packages.
#' @param ... additional keyword arguments indicating the rows corresponding to
#' elements in the seurat object metadata that should be filtered by.
#' For example, seurat_clusters = c(1, 9, 10) will filter the cells to only
#' those in seurat clusters 1, 9, and 10. Or, if there is some column named
#' "sample", then sample = c("A", "B") will conduct a run only on those
#' corresponding cells. A useful application for this is to run the function
#' many times, once for each sample, and plot them all together with [APOTCPlot]
#' down the line.
#' @param extra_filter character. An additional string that should be formatted
#' *exactly* like a statement one would pass into [dplyr::filter] that does
#' *additional* filtering to cells in the seurat object - on top of the other
#' keyword arguments - based on the metadata. This means that it will be
#' logically AND'ed with any keyword argument filters. This is a more flexible
#' alternative / addition to the filtering keyword arguments. For example, if
#' one wanted to filter by the length of the amino acid sequence of TCRs, one
#' could pass in something like `extra_filter = "nchar(CTaa) - 1 > 10"`. When
#' involving characters, ensure to enclose with single quotes.
#' @param run_id character. This will be the ID associated with the data of a
#' run, and will be used by other important functions like [APOTCPlot] and
#' [AdjustAPOTC]. Defaults to `NULL`, in which case the ID will be generated
#' in the following format:
#' `reduction_base;clonecall;keyword_arguments;extra_filter` where if
#' keyword arguments and extra_filter are underscore characters if there was
#' no input for the `...` and `extra_filter` parameters.
#' @param clone_scale_factor Dictates how much to scale each circle(between 0,1)
#' radius when converting from clonotype counts into circles that represent
#' individual clonotypes. The argument defaults to the character `"auto"`, and
#' if so, the most visually pleasing factor will be estimated
#' @param rad_scale_factor numeric between 0 and 1. This value decreases the
#' radius of the smallest clones by this scale factor. And the absolute value
#' of this decrease will be applied to all packed circles, effectively shrinking
#' all circles on the spot, and introduce more constant spacing in between
#' @param order_clones logical. Decides if the largest clone circles should be
#' near cluster centroids. This is highly recommended to be set to TRUE.
#' @param scramble_clones logical. Decides if the clone circles within each
#' cluster should be randomly scrambled when plotted. Note that solely this
#' argument or order_clones may be `TRUE` at once.
#' @param try_place If `TRUE`, always minimizes distance from a newly placed
#' circle to the origin in the circle packing algorithm
#' @param repulse If `TRUE`, will attempt to push overlapping clusters away from
#' each other
#' @param repulsion_threshold numeric. The radius that clonal circle clusters
#' overlap is acceptable when repulsing
#' @param repulsion_strength numeric. The smaller the value the less the
#' clusters repulse each other per iteration, and vice versa
#' @param max_repulsion_iter integer. The number of repulsion iterations
#' @param override logical. If `TRUE`, will override any existing
#' APackOfTheClones run data with the same parameters quietly
#' @param verbose logical. Decides if visual cues are displayed to the R console
#' of the progress
#'
#' @details
#' All APackOfTheClones run data is stored in the Seurat object under
#' `seurat_object@misc$APackOfTheClones`, which is a list of S4 objects of the
#' type "ApotcData", with each element corresponding to a unique run. The id of
#' each run is the name of each element in the list. The user
#' ***really shouldn't*** manually modify anything in the list as it may cause
#' unexpected behavior with many other functions.
#'
#' Additionally, it also logs a seurat command associated with the run in the
#' `seurat_object@commands` slot as a "SeuratCommand" object (from Seurat),
#' where the name of the object in the list is formatted as `RunAPOTC.run_id`.
#' Again, the user should not modify anything in these objects as they are used
#' by some related functions, mainly [AdjustAPOTC].
#'
#' @return A modified version of the input seurat object, which harbors data
#' necessary for visualizing the clonal expansion of the cells with [APOTCPlot]
#' and has a friendly user interface to modify certain attributes with
#' [AdjustAPOTC].
#'
#' @seealso [APOTCPlot], [AdjustAPOTC], [getApotcDataIds]
#'
#' @export
#'
#' @examples
#' data("combined_pbmc")
#'
# ' # this is the recommended approach to use a custom run_id with default params
# ' combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "default")
# '
# ' # here's a seperate run with some filters to the meta data, where
# ' # `orig.ident` is a custom column in the example data. Notice that it is not
# ' # a `RunAPOTC` parameter but a user keyword argument
# ' combined_pbmc <- RunAPOTC(
# '     combined_pbmc, run_id = "sample17", orig.ident = c("P17B", "P17L")
# ' )
# '
# ' # the exact same thing can be achieved with the `extra_filter` parameter
# ' combined_pbmc <- RunAPOTC(
# '     combined_pbmc,
# '     run_id = "sample17",
# '     extra_filter = "substr(orig.ident, 2, 3) == '17'",
# '     override = TRUE
# ' )
#'
RunAPOTC <- function(
    seurat_obj,
    reduction_base = "umap",
    clonecall = "strict",
    ...,
    extra_filter = NULL,
    run_id = NULL,

    clone_scale_factor = "auto",
    rad_scale_factor = 0.95,
    order_clones = TRUE,
    scramble_clones = FALSE,
    try_place = FALSE,
    repulse = TRUE,
    repulsion_threshold = 1,
    repulsion_strength = 1,
    max_repulsion_iter = 20L,

    override = TRUE,
    verbose = TRUE
) {
    call_time <- Sys.time()
    if (verbose) message("Initializing APOTC run...")

    # compute inputs
    reduction_base <- attempt_correction(reduction_base)
    clonecall <- .theCall(seurat_obj@meta.data, clonecall)
    varargs_list <- list(...)

    if (should_estimate(clone_scale_factor)) {
        clone_scale_factor <- estimate_clone_scale_factor(seurat_obj, clonecall)
        if (verbose) message(paste(
            "* Setting `clone_scale_factor` to",
            round(clone_scale_factor, digits = 3)
        ))
    }

    metadata_filter_string <- parse_to_metadata_filter_str(
        metadata_filter = extra_filter, varargs_list = varargs_list
    )

    obj_id <- infer_object_id_if_needed(
        args = hash::hash(as.list(environment())), varargs_list = varargs_list
    )

    RunAPOTC_parameter_checker(hash::hash(as.list(environment())))

    if (verbose) {
        message(paste("* id for this run:", obj_id))
        if (override && containsApotcRun(seurat_obj, obj_id))
            message("* overriding results of the previous run")
    }

    # run the packing algos
    apotc_obj <- ApotcData(
        seurat_obj,
        metadata_filter_string,
        clonecall,
        reduction_base,
        clone_scale_factor,
        rad_scale_factor
    )

    if (verbose) message("\nPacking clones into clusters")

    apotc_obj <- circlepackClones(
        apotc_obj,
        ORDER = order_clones,
        scramble = scramble_clones,
        try_place = try_place,
        verbose = verbose
    )

    if (repulse) {
        apotc_obj <- repulseClusters(
            apotc_obj,
            repulsion_threshold,
            repulsion_strength,
            max_repulsion_iter,
            verbose
        )
    }

    # store the apotc object in the correct slot with the correct id
    seurat_obj <- setApotcData(seurat_obj, obj_id, apotc_obj)

    # TODO this could be simplified, something with looking at the grandparent
    # environment frame with n = 2
    command_obj <- make_apotc_command(call_time)
    seurat_obj <- log_seurat_command(
        seurat_obj = seurat_obj,
        command_obj = command_obj,
        id = obj_id
    )

    if (verbose) print_completion_time(call_time, newline = TRUE)
    seurat_obj
}

RunAPOTC_parameter_checker <- function(args) {

    if (!is_seurat_object(args[["seurat_obj"]])) {
        stop("`seurat_obj` must be a Seurat object.")
    }

    # Check if reduction_base is character of length 1
    if (!is.character(args[["reduction_base"]]) || length(args[["reduction_base"]]) != 1) {
        stop("`reduction_base` must be a character of length 1.")
    }

    # Check if clonecall is character of length 1
    if (!is.character(args[["clonecall"]]) || length(args[["clonecall"]]) != 1) {
        stop("`clonecall` must be a character of length 1.")
    }

    # Check if extra_filter is character or NULL
    if (!is.null(args[["extra_filter"]]) && (!is.character(args[["extra_filter"]]) || length(args[["extra_filter"]]) != 1)) {
        stop("`extra_filter` must be a character or NULL of length 1.")
    }

    # Check if run_id is character or numeric or NULL
    if (!is.null(args[["run_id"]])) {
        if ((!is.character(args[["run_id"]]) && !is.numeric(args[["run_id"]])) || length(args[["run_id"]]) != 1) {
            stop("`run_id` must be a character or numeric of length 1.")
        }
    }

    # Check if clone_scale_factor is numeric of length 1
    if (!is.numeric(args[["clone_scale_factor"]]) || length(args[["clone_scale_factor"]]) != 1) {
        stop("`clone_scale_factor` must be a numeric value of length 1.")
    }

    # Check if rad_scale_factor is numeric of length 1
    if (!is.numeric(args[["rad_scale_factor"]]) || length(args[["rad_scale_factor"]]) != 1) {
        stop("`rad_scale_factor` must be a numeric value of length 1.")
    }

    # Check if order_clones and scramble_clones are logical
    if (!is.logical(args[["order_clones"]]) || length(args[["order_clones"]]) != 1 ||
        !is.logical(args[["scramble_clones"]]) || length(args[["scramble_clones"]]) != 1) {
        stop("`order_clones` and `scramble_clones` must be logical values of length 1.")
    }

    # Check if try_place is logical of length 1
    if (!is.logical(args[["try_place"]]) || length(args[["try_place"]]) != 1) {
        stop("`try_place` must be a logical value of length 1.")
    }

    # Check if repulse is logical of length 1
    if (!is.logical(args[["repulse"]]) || length(args[["repulse"]]) != 1) {
        stop("`repulse` must be a logical value of length 1.")
    }

    # Check if repulsion_threshold is numeric of length 1
    if (!is.numeric(args[["repulsion_threshold"]]) || length(args[["repulsion_threshold"]]) != 1) {
        stop("`repulsion_threshold` must be a numeric value of length 1.")
    }

    # Check if repulsion_strength is numeric of length 1
    if (!is.numeric(args[["repulsion_strength"]]) || length(args[["repulsion_strength"]]) != 1) {
        stop("`repulsion_strength` must be a numeric value of length 1.")
    }

    # Check if max_repulsion_iter is an integer of length 1
    if (!is.integer(args[["max_repulsion_iter"]]) || length(args[["max_repulsion_iter"]]) != 1) {
        stop("`max_repulsion_iter` must be an integer value of length 1.")
    }

    # Check if override is logical of length 1
    if (!is.logical(args[["override"]]) || length(args[["override"]]) != 1) {
        stop("`override` must be a logical value of length 1.")
    }

    # Check if verbose is logical of length 1
    if (!is.logical(args[["verbose"]]) || length(args[["verbose"]]) != 1) {
        stop("`verbose` must be a logical value of length 1.")
    }

    # regular tests

	if (args[["clone_scale_factor"]] <= 0 || args[["clone_scale_factor"]] > 1) {
		stop("`clone_scale_factor` has to be a positive real number in (0, 1]")
	}

    if (args[["rad_scale_factor"]] <= 0 || args[["rad_scale_factor"]] > 1) {
		stop("`rad_scale_factor` has to be a positive real number in (0, 1]")
	}

	if (args[["order_clones"]] == args[["scramble_clones"]]) {
		stop(paste(
            "`order_clones` and `scramble_clones` cannot both be",
            args[["order_clones"]]
        ))
	}

    if (args[["repulse"]]) {
        if (args[["repulsion_threshold"]] <= 0) {
            stop("`repulsion_threshold` has to be a positive number")
        }
        if (args[["repulsion_strength"]] <= 0) {
            stop("`repulsion_strength` has to be a positive number")
        }
        if (args[["max_repulsion_iter"]] <= 0) {
            stop("`max_repulsion_iter` has to be a positive number")
        }
    }

    if (!args[["override"]] && containsApotcRun(args[["seurat_obj"]], args[["obj_id"]])) {
        stop(paste(
            "An APackOfTheClones run with the the parameters", args[["obj_id"]],
            "appears to already have been ran. If this is a mistake,",
            "set the `override` argument to `TRUE` and re-run."
        ))
    }

	# TODO more checks of the filtering conditions and more tests

    check_filtering_conditions(args)
}

# assumes varargs_list is present
check_filtering_conditions <- function(args) {
    if (is_empty(args$varargs_list)) return()
    metadata_cols <- names(args$seurat_obj@meta.data)
    all_formals <- get_processed_argnames(3)
    for (argname in names(args$varargs_list)) {
        if (argname %in% metadata_cols) next
        stop(paste(
            "colname:", argname,
            "not found in the seurat object metadata.",
            "did you mean this to be a subsetting named argument?",
            "if not, did you mean to use the argument:",
            closest_word(argname, all_formals)
        ))
    }
}