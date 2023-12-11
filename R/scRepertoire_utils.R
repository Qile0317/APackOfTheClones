.convertClonecall <- function(x) {

  clonecall_dictionary <- hash::hash(
    "gene" = "CTgene",
		"genes" = "CTgene",
		"ctgene" = "CTgene",
		"ctstrict" = "CTstrict",
		"nt" = "CTnt",
		"nucleotide" = "CTnt",
		"nucleotides" = "CTnt",
		"ctnt" = "CTnt",
		"aa" = "CTaa",
		"amino" = "CTaa",
		"ctaa" = "CTaa",
		"gene+nt" = "CTstrict",
		"strict" = "CTstrict",
		"ctstrict" = "CTstrict"
	)

	x <- tolower(x)

	if (!is.null(clonecall_dictionary[[x]])) {
		return(clonecall_dictionary[[x]])
	}

	stop(paste(
		"invalid input cloneCall, did you mean: '",
		closest_word(
			x,
			c(names(clonecall_dictionary),
			  unname(hash::values(clonecall_dictionary)))
		),
		"'?",
		sep = ""
	))
}
