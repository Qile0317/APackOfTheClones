#' @title
#' Artificially generated Seurat object
#'
#' @description
#' A generated 'SeuratObject' of a small single-sample sc-RNAseq experiment.
#' Has a corresponding T-cell receptor library generated from
#' single cell immune profiling, named `"mini_clonotype_data"`
#'
#' @usage data("mini_seurat_obj")
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
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: PCA, UMAP, and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#'
#' @seealso [mini_clonotype_data()]
"mini_seurat_obj"

#' @title
#' Artificially generated T cell receptor library
#' 
#' @keywords internal
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' A generated dataframe of a T-cell receptor (TCR) library generated from
#' single cell immune profiling. It is a subset the full dataframe
#' which would usually have up to 18 columns containing different data,
#' because the intended purpose of this object is to test various functions
#' in 'APackOfTheClones'. The dataframe compliments `mini_seurat_obj` and
#' can be integrated into it with `integrate_tcr`.
#'
#' @usage data("mini_clonotype_data")
#'
#' @format `data.frame`
#' A data frame with 80 rows and 2 columns:
#' \describe{
#'   \item{barcode}{barcodes corresponding to each sequenced cell}
#'   \item{raw_clonotype_id}{clonotype information for each cell}
#' }
#'
#' @details Note that the clonotypes in the `raw_clonotype_id` column
#' actually do not contain all of clonotype`1`...clonotype`n`
#'
#' @seealso [mini_seurat_obj()]
"mini_clonotype_data"

#' @title
#' Multi-sampled T-cell seurat object with integrated TCR library
#'
#' @description
#' Generated with [combineSeuratExpression], more specifically, with:
#'
#' `r get(data(combined_pbmc))@commands[["combineSeuratExpression"]]@call.string`
#'
#' @usage data("combined_pbmc")
#'
#' @format A seurat object
#'
"combined_pbmc"
