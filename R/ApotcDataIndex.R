# script to manage the interface for accessing the apotc data
# could also put it in ApotcData.R

.defaultApotcDataObjId <- "__all__"
utils::globalVariables(c(".defaultApotcDataSample"))

# from the input of RunAPOTC, convert the condition to a call to be put in
# dpylr::filter for the metadata
parse_to_metadata_filter <- function(filter_samples, filter_ID, metadata_filter) {
    if (all(is.null(c(filter_samples, filter_ID, metadata_filter)))) {
        return(NULL)
    }
    # TODO
}

# from the input of RunAPOTC, convert the condition to the apotc data sample id where
# its stored under @misc[["APackOfTheClones"]][[id]]
parse_to_object_id <- function(reduction_base, clonecall, ...) {
	if (is.null(filter_string)) {
		return(.defaultApotcDataSample)
	}
	# TODO
}

# make user getters

# new format, there will be a list of apotc objects in the seurat@misc slot. the list will be named apotc.
# each is dependent on reduction/samples and within the list there will be named elements for each reduction/sample combo
# and make it optional in RunAPOTC if this should be stored. APOTCPlot will then be able to have the apotc obj slot input
# alternatively the sample/id configuration.
