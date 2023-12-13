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
#' dimensional  reduction (`reduction_base`) coordinates.
#'
#' The parameter `extra_filter` along with an unlimited number of additional
#' keyword arguments can be used to filter the cells by certain conditions
#' in the metadata, and new results will be stored in addition to other runs
#' the users may have done.
#'
#' If the user wishes to manually customize/fix the expansion plot
#' generated, the circular packing information can be modified
#' with the [AdjustAPOTC] function.
#'
#' @param seurat_obj Seurat object with one or more dimension reductions. Must
#' already have been integrated with a T cell library via
#' `integrate_tcr(seurat_obj, tcr_df)`
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
#' keyword arguments - based on the metadata. This is a more flexible
#' alternative / addition to the filtering keyword arguments. For example, if
#' one wanted to filter by the length of the amino acid sequence of TCRs, one
#' could pass in something like `extra_filter = "nchar(CTaa) - 1 > 10"`
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
#' Each APackOfTheClones run is uniquely identified by the parameters
#' `reduction_base`, `clonecall`, `extra_filter`, and any additional keywords
#' passed to filter the metadata. Each distinct run is stored in the seurat
#' object under the `@misc` slot in a list named `"APackOfTheClones"`, where
#' each list element is an instance of the [ApotcData] class. Each element is
#' named with an object id which is generated from the aforementioned parameters
#' stored. The user can but is recommended heavily to not modify anything under
#' the @misc$APackOfTheClones manually as it may cause unexpected behavior.
#'
#' @return A modified version of the input seurat object, which harbors data
#' necessary for visualizing the clonal expansion of the cells with [APOTCPlot]
#' and has a friendly user interface to modify certain attributes with
#' [AdjustAPOTC].
#'
#' @seealso [APOTCPlot], [AdjustAPOTC]
#'
#' @export
#'
#' @examples
#' data("combined_pbmc")
#' combined_pbmc <- RunAPOTC(combined_pbmc, verbose = FALSE)
#'
RunAPOTC <- function(
    seurat_obj,
    reduction_base = "umap",
    clonecall = "strict",
    ...,
    extra_filter = NULL,

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

    if (should_estimate(clone_scale_factor)) {
        clone_scale_factor <- estimate_clone_scale_factor(seurat_obj, clonecall)
        if (verbose) message(paste(
            "* Setting `clone_scale_factor` to",
            round(clone_scale_factor, digits = 3)
        ))
    }

    #clonecall <- scRepertoire:::.theCall(clonecall, seurat_obj@meta.data)
    clonecall <- .convertClonecall(clonecall)

    metadata_filter_string <- parse_to_metadata_filter_str(
        metadata_filter = extra_filter, varargs_list = list(...)
    )

    obj_id <- parse_to_object_id(
        reduction_base = reduction_base, clonecall =  clonecall,
        varargs_list = list(...), metadata_filter = extra_filter
    )

    RunAPOTC_parameter_checker(hash::hash(as.list(environment())))

    if (verbose) message(paste("* id for this run:", obj_id, "\n"))

    # run the packing algos
    apotc_obj <- ApotcData(
        seurat_obj, metadata_filter_string, clonecall, reduction_base,
        clone_scale_factor, rad_scale_factor
    )

    if (verbose) message("Packing clones into clusters")

    apotc_obj <- circlepackClones(
        apotc_obj, order_clones, scramble_clones, try_place, verbose
    )

    if (repulse) {
        apotc_obj <- repulseClusters(
            apotc_obj, repulsion_threshold, repulsion_strength,
            max_repulsion_iter, verbose
        )
    }

    # store the apotc object in the correct slot with the correct id
    seurat_obj <- setApotcData(seurat_obj, obj_id, apotc_obj)

    seurat_obj <- log_and_index_command(
        seurat_obj, "RunAPOTC", command_obj = make_apotc_command(call_time)
    )

    if (verbose) print_completion_time(call_time, newline = TRUE)
    seurat_obj
}

RunAPOTC_parameter_checker <- function(args) {

	if (is.null(args[["seurat_obj"]]@reductions[[attempt_correction(args[["reduction_base"]])]])) {
		stop(paste(
			"No", args[["reduction_base"]], "reduction found on the seurat object,",
			"ensure the the reduction has been computed. Otherwise, did you",
			"mean:", closest_word(args[["reduction_base"]], c("umap", "tsne", "pca"))
		))
	}

	if (args[["clone_scale_factor"]] <= 0 || args[["clone_scale_factor"]] > 1) {
		stop("`clone_scale_factor` has to be a positive number in (0, 1]")
	}

    if (args[["rad_scale_factor"]] <= 0 || args[["rad_scale_factor"]] > 1) {
		stop("`rad_scale_factor` has to be a positive number in (0, 1]")
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

    if (args[["override"]] && containsApotcRun(args[["seurat_obj"]], args[["obj_id"]])) {
        stop(paste(
            "An APackOfTheClones run with the the parameters", args[["obj_id"]],
            "appears to already have been ran. If this is a mistake,",
            "set the `override` argument to `FALSE` and re-run."
        ))
    }

	# TODO more checks of the filtering conditions
}
