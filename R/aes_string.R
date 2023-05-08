library(ggplot2)
library(rlang)

# the following code is copied and modified from ggplot2
# under the R/aes.R file (under MIT license, like this package)
rename_aes <- function(x) {
  names(x) <- ggplot2::standardise_aes_names(names(x))
  duplicated_names <- names(x)[duplicated(names(x))]
  if (length(duplicated_names) > 0L) {
    return(NULL)
  }
  x
}

aes_string <- function(x, y, ...) {
  mapping <- list(...)
  if (!missing(x)) mapping["x"] <- list(x)
  if (!missing(y)) mapping["y"] <- list(y)
  
  mapping <- lapply(mapping, function(x) {
    if (is.character(x)) {
      x <- rlang::parse_expr(x)
    }
  })
  structure(rename_aes(mapping), class = "uneval")
}