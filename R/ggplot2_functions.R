# the following code is copied and modified from ggplot2 under the MIT license

apotc_rename_aes <- function(x) {
  names(x) <- ggplot2::standardise_aes_names(names(x))
  duplicated_names <- names(x)[duplicated(names(x))]
  if (length(duplicated_names) > 0L) {
    return(NULL)
  }
  x
}

apotc_aes_string <- function(x, y, ...) {
  mapping <- list(...)
  if (!missing(x)) mapping["x"] <- list(x)
  if (!missing(y)) mapping["y"] <- list(y)

  mapping <- lapply(mapping, function(x) {
    if (is.character(x)) {
      x <- rlang::parse_expr(x)
    }
  })
  structure(apotc_rename_aes(mapping), class = "uneval")
}
