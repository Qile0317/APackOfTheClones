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
#' `scRepertoire::combineExpression`.
#' @param reduction_base character. The seurat reduction to base the clonal
#' expansion plotting on. Defaults to `'umap'` but can be any reduction present
#' within the reductions slot of the input seurat object, including custom ones.
#' If `'pca'``, the cluster coordinates will be based on PC1 and PC2.
#' However, generally APackOfTheClones is used for displaying UMAP and
#' occasionally t-SNE versions to intuitively highlight clonal expansion.
#' @param clonecall character. The column name in the seurat object metadata to
#' use. See `scRepertoire` documentation for more information about this
#' parameter that is central to both packages.
#' @param ... additional "subsetting" keyword arguments indicating the rows
#' corresponding to elements in the seurat object metadata that should be
#' filtered by. E.g., `seurat_clusters = c(1, 9, 10)` will filter the cells to
#' those in the `seurat_clusters` column with any of the values 1, 9, and 10.
#' Unfortunately, column names in the seurat object metadata cannot
#' conflict with the keyword arguments. ***MAJOR NOTE*** if any subsetting
#' keyword arguments are a *prefix* of any preceding argument names (e.g. a
#' column named `reduction` is a prefix of the `reduction_base` argument)
#' R will interpret it as the same argument unless *both* arguments
#' are named. Additionally, this means any subsequent arguments *must* be named.
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
#'
#' `reduction_base;clonecall;keyword_arguments;extra_filter`
#'
#' where if keyword arguments and extra_filter are underscore characters if
#' there was no input for the `...` and `extra_filter` parameters.
#' @param clone_scale_factor Dictates how much to scale each circle(between 0,1)
#' radius when converting from clonotype counts into circles that represent
#' individual clonotypes. The argument defaults to the character `"auto"`, and
#' if so, the most visually pleasing factor will be estimated.
#' @param rad_scale_factor numeric between 0 and 1. This value decreases the
#' radius of the smallest clones by this scale factor. And the absolute value
#' of this decrease will be applied to all packed circles, effectively shrinking
#' all circles on the spot, and introduce more constant spacing in between.
#' @param order_clones logical. Decides if the largest clone circles should be
#' near cluster centroids. This is highly recommended to be set to TRUE for
#' increased intuitiveness of the visualization, as resulting plots tend to
#' give an improved impression of the proportion of expanded clones. If
#' `FALSE,` will randomly scramble the positions of each circle. For the sake
#' of being replicable, a random seed is recommended to be set with [set.seed].
#' @param try_place If `TRUE`, always minimizes distance from a newly placed
#' circle to the origin in the circle packing algorithm.
#' @param repulse If `TRUE`, will attempt to push overlapping clusters away from
#' each other.
#' @param repulsion_threshold numeric. The radius that clonal circle clusters
#' overlap is acceptable when repulsing.
#' @param repulsion_strength numeric. The smaller the value the less the
#' clusters repulse each other per iteration, and vice versa.
#' @param max_repulsion_iter integer. The number of repulsion iterations.
#' @param override logical. If `TRUE`, will override any existing
#' APackOfTheClones run data with the same `run_id`.
#' @param verbose logical. Decides if visual cues are displayed to the R console
#' of the progress.
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
#' # this is the recommended approach to use a custom run_id with default params
#' combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "default", verbose = FALSE)
#'
#' # here's a seperate run with some filters to the meta data, where
#' # `orig.ident` is a custom column in the example data. Notice that it is not
#' # a `RunAPOTC` parameter but a user keyword argument
#' combined_pbmc <- RunAPOTC(
#'     combined_pbmc, run_id = "sample17", orig.ident = c("P17B", "P17L"),
#'     verbose = FALSE
#' )
#'
#' # the exact same thing can be achieved with the `extra_filter` parameter
#' combined_pbmc <- RunAPOTC(
#'     combined_pbmc,
#'     run_id = "sample17",
#'     extra_filter = "substr(orig.ident, 2, 3) == '17'",
#'     override = TRUE,
#'     verbose = FALSE
#' )
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
    try_place = FALSE,

    repulse = TRUE,
    repulsion_threshold = 1,
    repulsion_strength = 1,
    max_repulsion_iter = 20L,

    override = FALSE,
    verbose = TRUE
) {
    # setup and check inputs
    call_time <- Sys.time()
    varargs_list <- list(...)
    RunAPOTC_partial_arg_checker(varargs_list)
    if (verbose) message("Initializing APOTC run...")

    # compute inputs
    reduction_base <- attempt_correction(seurat_obj, reduction_base)
    clonecall <- .theCall(seurat_obj@meta.data, clonecall)
    
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
        if (override && containsApotcRun(seurat_obj, obj_id)) {
            message("* overriding results of the previous run")
        }
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

# # a super simple regression of manually determined clone scale based on cell
# # count a much more complicated model can use other facts about the seurat
# # obj to improve how visually plesant is it. Some overlap factor could 
# # probably be eestimated with the raw clone counts
# cell_count <- c(80, 365, 2500)
# desirable_factor <- c(1, 0.3, 0.2)
# plot(cell_count, desirable_factor)

estimate_clone_scale_factor <- function(seurat_obj, clonecall) {
	num_clones <- count_clones(seurat_obj, clonecall)

	if (num_clones <= 365) {
		approx_clone_scale_factor <- (-0.002456 * num_clones) + 1.196491
	} else {
		approx_clone_scale_factor <- (-4.684e-05 * num_clones) + 3.171e-01
	}
	
	bound_num(approx_clone_scale_factor, lowerbound = 0.05, upperbound = 1)
}

RunAPOTC_partial_arg_checker <- function(varargs_list = list()) {

    args <- get_parent_func_args()

    check_apotc_identifiers(args)
    check_filtering_conditions(args)

    if (!should_estimate(args$clone_scale_factor)) {
        typecheck(args$clone_scale_factor, is_a_positive_numeric)
    }

    typecheck(args$rad_scale_factor, is_a_positive_numeric)

    typecheck(args$try_place, is_a_logical)

    check_repulsion_params(args)

    typecheck(args$override, is_a_logical)

    typecheck(args$verbose, is_a_logical)

}

check_apotc_identifiers <- function(args) {

    if (!is_seurat_object(args[["seurat_obj"]])) {
        stop(call. = FALSE, "`seurat_obj` must be a Seurat object.")
    }

    # TODO nulls shouldnt be allowed in some cases
    typecheck(args$reduction_base, is_a_character, is.null)
    typecheck(args$clonecall, is_a_character, is.null)
    typecheck(args$extra_filter, is_character, is.null)

}

check_filtering_conditions <- function(args, frame_level = 2) {

    if (is_empty(args$varargs_list)) return()

    metadata_cols <- names(args$seurat_obj@meta.data)
    all_formals <- get_processed_argnames(frame_level)

    for (argname in names(args$varargs_list)) {
        if (argname %in% metadata_cols) next
        stop(call. = FALSE,
            "should `", argname, "` be a function argument? ",
            "If so, did you mean `", closest_word(argname, all_formals), "`? ",
            "Otherwise, should `", argname, "` be a subsetting argument? ",
            "If so, did you mean `", closest_word(argname, metadata_cols), "`?"
        )
    }
}

check_repulsion_params <- function(args) {
    typecheck(args$repulse, is_a_logical)
    if (!args$repulse) return()
    typecheck(args$repulsion_threshold, is_a_positive_numeric)
    typecheck(args$repulsion_strength, is_a_positive_numeric)
    typecheck(args$max_repulsion_iter, is_a_positive_integer)
}

RunAPOTC_parameter_checker <- function(args) {

    if (!args$override && containsApotcRun(args$seurat_obj, args$obj_id)) {
        stop(call. = FALSE, paste(
            "An APackOfTheClones run with the the parameters", args$obj_id,
            "appears to already have been ran. If this is intended,",
            "set the `override` argument to `TRUE` and re-run."
        ))
    }

	# TODO more checks
}
