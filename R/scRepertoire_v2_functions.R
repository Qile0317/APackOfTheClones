# A variety of scRepertoire v2 functions copied during beta testing,
# all permissions granted by the owner + both authors

.theCall <- function(df, x, check.df = TRUE) {
    x <- .convertClonecall(x)
    if(check.df) {
      if(inherits(df, "list") & !any(colnames(df[[1]]) %in% x)) {
        stop("Check the clonal variabe (cloneCall) being used in the function, it does not appear in the data provided.")
      } else if (inherits(df, "data.frame") & !any(colnames(df) %in% x)) {
        stop("Check the clonal variabe (cloneCall) being used in the function, it does not appear in the data provided.")
      }
    }
    return(x)
}

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