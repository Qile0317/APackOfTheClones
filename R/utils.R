progress_bar <- function (x, max = 100) {
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
