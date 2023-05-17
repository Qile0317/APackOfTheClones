progress_bar <- function (x = 0, max = 100) {
  percent <- 100 * (x / max)
  cat(sprintf(
    '\r[%-50s] %d%%',
    paste(rep('=', percent * 0.5), collapse = ''),
    floor(percent))
  )
}

get_num_clusters <- function(seurat_obj) {
  length(levels(seurat_obj@meta.data[["seurat_clusters"]]))
}

isnt_empty <- function(inp) {
  !identical(inp, list())
}

isnt_na <- function(inp) {
  !any(is.na(inp))
}

isnt_empty_nor_na <- function(inp) {
  isnt_empty(inp) && isnt_na(inp)
}

#load_packages <- function(...) {
#  packages <- list(...)
#  n <- length(packages)
#  for (i in 1:n) {
#    suppressPackageStartupMessages(require(
#      packages[[i]], character.only = TRUE
#    ))
#    progress_bar(i, n)
#  }
#}
##removed from stackoverflow advice
