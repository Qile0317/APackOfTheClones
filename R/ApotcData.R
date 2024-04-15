methods::setClass(
	Class = "ApotcData",
	slots = c(
		reduction_base = 'character',
		clonecall = 'character',
		metadata_filter_string = 'character',
		# could add ident_col but not really needed?

		clusters = 'list',
		centroids = 'list',
		clone_sizes = 'list',
		idents = "character",

		clone_scale_factor = 'numeric',
		rad_scale_factor = 'numeric',
		cluster_colors = 'character',
		labels = 'character',
		label_coords = 'list'
	)
)

# initialize the ApotcData object with all slots except clusters done
# TODO a subsetted seurat object will probably have wrong colors
ApotcData <- function(
	seurat_obj,
	alt_ident,
	metadata_filter_condition,
	clonecall,
	reduction_base,
	clone_scale_factor,
	rad_scale_factor
) {

	apotc_obj <- initializeApotcData(
		seurat_obj, clonecall, metadata_filter_condition, reduction_base,
		clone_scale_factor, rad_scale_factor
	)
	
	seurat_obj <- set_meta_ident_col(seurat_obj, alt_ident)

	apotc_obj <- addIdentsLabelsColors(apotc_obj, seurat_obj)

	if (!identical(metadata_filter_condition, "")) {
		seurat_obj <- subsetSeuratMetaData(
			seurat_obj, metadata_filter_condition
		)
	}

	seurat_obj <- ident_into_seurat_clusters(seurat_obj)

	apotc_obj %>%
		addRawClusteredCloneSizes(seurat_obj) %>%
		addInitialCentroids(seurat_obj)
}

initializeApotcData <- function(
	seurat_obj, clonecall, filter_string, reduction_base,
	clone_scale_factor, rad_scale_factor
) {
	methods::new(
		Class = "ApotcData",

		reduction_base = reduction_base,
		clonecall = clonecall,
		metadata_filter_string = filter_string,

		clusters = list(),
		centroids = list(),
		clone_sizes = list(),
		idents = character(0),

		clone_scale_factor = clone_scale_factor,
		rad_scale_factor = rad_scale_factor,
		cluster_colors = character(0),
		labels = character(0),
		label_coords = list()
	)
}

# terrible hack fix to set seurat_clusters column to active ident or alt_ident
# so that the cluster counting is consistent - assumes alt_ident is valid
# the output is a seurat obj with a new column in the metadata called
# __active.ident__ that is guaranteed to be a factor with unguaranteed level type
set_meta_ident_col <- function(seurat_obj, alt_ident) {

    if (is.null(alt_ident)) {
        seurat_obj@meta.data[["__active.ident__"]] <- seurat_obj@active.ident
        return(seurat_obj)
    }

    seurat_obj@meta.data[["__active.ident__"]] <- as.factor(
        seurat_obj@meta.data[[alt_ident]]
    )

    seurat_obj
}

# finish the idents and num_clusters slot
# apotc_obj must be a product of initializeApotcData
# seurat_obj must be a product of set_meta-ident_col
addIdentsLabelsColors <- function(apotc_obj, seurat_obj) {

	apotc_obj <- set_idents(
		apotc_obj, levels(seurat_obj@meta.data[["__active.ident__"]])
	)

	apotc_obj <- set_cluster_colors(
		apotc_obj, gg_color_hue(get_num_clusters(apotc_obj))
	)

	if (have_default_idents(apotc_obj)) {
		return(set_labels(apotc_obj, gen_labels(get_num_clusters(apotc_obj))))
	}

	set_labels(apotc_obj, get_idents(apotc_obj))
}

# important function to be ran after setting meta ident col
# to merge that temporary column into seurat_clusters column
ident_into_seurat_clusters <- function(seurat_obj) {
    idents <- seurat_obj@meta.data[["__active.ident__"]]
    seurat_obj@meta.data[["__active.ident__"]] <- NULL
    seurat_obj@meta.data$seurat_clusters <- idents
    seurat_obj
}

# seurat obj has seurat clusters column with correct idents
addRawClusteredCloneSizes <- function(apotc_obj, seurat_obj) {
	set_raw_clone_sizes(apotc_obj, count_raw_clone_sizes(
		seurat_obj, get_idents(apotc_obj), get_clonecall(apotc_obj)
	))
}

# assume apotc_obj has correct ident levels
addInitialCentroids <- function(apotc_obj, seurat_obj) {
	initial_centroids <- get_cluster_centroids(
		seurat_obj, get_reduction_base(apotc_obj), get_idents(apotc_obj)
	)

	apotc_obj@centroids <- apotc_obj@label_coords <- initial_centroids
	apotc_obj
}

# pack the clones assuming centroids are present
circlepackClones <- function(apotc_obj, ORDER, try_place, verbose) {

	apotc_obj@clusters <- pack_into_clusterlists(
		sizes = get_processed_clone_sizes(apotc_obj),
		centroids = get_centroids(apotc_obj),
		num_clusters = get_num_clusters(apotc_obj),
		rad_decrease = get_rad_decrease(apotc_obj),
		ORDER = ORDER,
		scramble = !ORDER,
		try_place = try_place,
		verbose = verbose
	)

	# see which elemens of sizes are empty and set corresponding elements empty
	for (i in seq_len(get_num_clusters(apotc_obj))) {
		if (isnt_empty(apotc_obj@clusters[[i]])) {
			next
		}
		apotc_obj@centroids[[i]] <- list()
		apotc_obj@label_coords[[i]] <- list()
	}

	apotc_obj
}

repulseClusters <- function(
	apotc_obj, repulsion_threshold, repulsion_strength, max_repulsion_iter,
	verbose
) {
	repulsed_clusters <- get_repulsed_clusterlists(
	    packed_clusters = get_clusterlists(apotc_obj),
	    repulsion_threshold = repulsion_threshold,
		repulsion_strength = repulsion_strength,
		max_repulsion_iter = max_repulsion_iter,
		verbose = verbose
	)

	setModifiedClusterlists(
		apotc_obj, modified_clusterlists = repulsed_clusters
	)
}

# function to modify the apotc_obj's relevant slots when modified clusterlists
# are introduced e.g. for cluster repulsion or relocation. This cannot be used
# for completely new irrelevant clusterlists, as the centroids and label_coords
# are modified correspondingly to the original clusters.
setModifiedClusterlists <- function(apotc_obj, modified_clusterlists) {

	modified_centroids <- read_centroids(modified_clusterlists)

	apotc_obj@label_coords <- move_coord_list_by_same_amount(
		coord_list = get_label_coords(apotc_obj),
		original_coord_list = get_centroids(apotc_obj),
		new_coord_list = modified_centroids
	)

	apotc_obj@clusters <- modified_clusterlists
	apotc_obj@centroids <- modified_centroids

	apotc_obj
}

convert_to_rad_decrease <- function(clone_scale_factor, rad_scale_factor) {
	clone_scale_factor * (1 - rad_scale_factor)
}

match_index <- function(apotc_obj, index) {

	varname <- deparse(substitute(index))

	if (is_integer(index)) {
		if(all(sapply(index, function(x) is_valid_cluster(apotc_obj, x)))) {
			return(index)
		}
		stop(call. = FALSE,
			"Some or all indices in `", varname, "` ", #FIXME likely wrong varname
			"are out of bounds of the APackOfTheClones Run."
		)
	}

	# assume character of names
	labels <- get_labels(apotc_obj)
	output <- integer(length(labels))

	for (i in seq_along(index)) {

		location <- which(labels == index[i]) # TODO check invalid labels

		if (length(location) == 0) {
			stop(call. = FALSE,
				"No label matched the input string"
			)
		}

		if (length(location) > 1) {
			warning(call. = FALSE,
				"* label '", index[i], "' ",
				"had multiple matches in the APackOfTheClones run, ",
				"using the first occurence at index ", location[1]
			)
		}

		output[i] <- location[1]
	}

	output
}

# checkers

is_valid_cluster <- function(apotc_obj, cluster_ind) {
	cluster_ind %>% is_bound_between(1, get_num_clusters(apotc_obj))
}

is_valid_nonempty_cluster <- function(apotc_obj, cluster_ind) {
	typecheck(cluster_ind, is_an_integer)
	is_valid_cluster(apotc_obj, cluster_ind) &&
		isnt_empty(get_clusterlists(apotc_obj)[[cluster_ind]])
}

have_default_idents <- function(apotc_obj) {
	identical(get_idents(apotc_obj), as.character(1:get_num_clusters(apotc_obj)))
}

# getters

get_reduction_base <- function(apotc_obj) {
	apotc_obj@reduction_base
}

get_clonecall <- function(apotc_obj) {
	apotc_obj@clonecall
}

get_metadata_filter_string <- function(apotc_obj) {
	apotc_obj@metadata_filter_string
}

get_idents <- function(apotc_obj) {
	apotc_obj@idents
}

get_ident_index <- function(apotc_obj, ident_lvl) {
	which(get_idents(apotc_obj) %in% ident_lvl)
}

get_clusterlists <- function(apotc_obj) {
	apotc_obj@clusters
}

get_centroids <- function(apotc_obj) {
	apotc_obj@centroids
}

get_raw_clone_sizes <- function(apotc_obj, as_hash = FALSE) {
	if (!as_hash) return(apotc_obj@clone_sizes)
	hash_from_tablelist(apotc_obj@clone_sizes)
}

# get all clone sizes as a single table
# TODO unsure what happens with no clones
get_aggregated_clone_sizes <- function(
	apotc_obj, sort_decreasing = NULL, get_top = NULL
) {
	clone_sizes <- aggregate_clone_sizes(
		get_raw_clone_sizes(apotc_obj), sort_decreasing
	)

	if (is.null(get_top)) return(clone_sizes)
	filter_top_clones(clone_sizes, get_top)
}

get_unique_clonotypes <- function(x) {
	unique(unlist(lapply(get_raw_clone_sizes(x), names)))
}

get_processed_clone_sizes <- function(apotc_obj) {
  raw_tabled_clone_sizes <- get_raw_clone_sizes(apotc_obj)
  processed_sizes <- init_list(get_num_clusters(apotc_obj), list())

  for (i in seq_len(get_num_clusters(apotc_obj))) {
    if (!is_empty_table(raw_tabled_clone_sizes[[i]])) {
      processed_sizes[[i]] <- apotc_obj@clone_scale_factor *
        sqrt(raw_tabled_clone_sizes[[i]])
    }
  }
  processed_sizes
}

get_num_clones <- function(apotc_obj) {
	sum(unlist(get_raw_clone_sizes(apotc_obj)))
}

get_num_clusters <- function(apotc_obj) length(get_idents(apotc_obj))

get_valid_num_clusters <- function(apotc_obj) {
	n <- 0
	for (cluster in apotc_obj@clusters) {
		if (isnt_empty(cluster)) {
			n <- n + 1
		}
	}
	n
}

get_clone_scale_factor <- function(apotc_obj) {
	apotc_obj@clone_scale_factor
}

get_rad_scale_factor <- function(apotc_obj) {
	apotc_obj@rad_scale_factor
}

get_rad_decrease <- function(apotc_obj) {
	convert_to_rad_decrease(
		clone_scale_factor = get_clone_scale_factor(apotc_obj),
		rad_scale_factor = get_rad_scale_factor(apotc_obj)
	)
}

get_cluster_colors <- function(apotc_obj) {
	apotc_obj@cluster_colors
}

get_labels <- function(apotc_obj) {
	apotc_obj@labels
}

get_label_coords <- function(apotc_obj) {
	apotc_obj@label_coords
}

# setters

set_metadata_filter_string <- function(apotc_obj, extra_filter) {
	apotc_obj@metadata_filter_string <- extra_filter
	apotc_obj
}

set_idents <- function(apotc_obj, x) {
	apotc_obj@idents <- as.character(x)
	apotc_obj
}

set_raw_clone_sizes <- function(apotc_obj, x) {
	apotc_obj@clone_sizes <- x
	apotc_obj
}

set_clone_scale_factor <- function(apotc_obj, x) {
	apotc_obj@clone_scale_factor <- x
	apotc_obj
}

set_clusterlists <- function(apotc_obj, x) {
	apotc_obj@clusters <- x
	apotc_obj
}

lapply_clusterlists <- function(apotc_obj, f) {
	set_clusterlists(apotc_obj, lapply(get_clusterlists(apotc_obj), f))
}

set_cluster_colors <- function(apotc_obj, x) {
	apotc_obj@cluster_colors <- x
	apotc_obj
}

set_labels <- function(apotc_obj, x) {
	apotc_obj@labels <- x
	apotc_obj
}
