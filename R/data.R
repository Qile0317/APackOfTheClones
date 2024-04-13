#' @title
#' Example Multi-sampled T-cell seurat object with integrated TCR library
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Generated with `scRepertoire::combineExpression`. To construct this object
#' from scratch, try:
#'
#' `r get(data(combined_pbmc))@commands$combineExpression@call.string`
#'
#' @usage data("combined_pbmc")
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata with a combined TCR contig list `from scRepertoire`}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: UMAP}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history, including the one used to create this object "combineExpression"}
#' }
#'
"combined_pbmc"
