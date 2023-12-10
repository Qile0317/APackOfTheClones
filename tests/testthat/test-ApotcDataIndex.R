test_that("parse_to_metadata_filter_str works", {
    # TODO
})

test_that("parse_to_object_id works for single varargs", {

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap", clonecall = "CTstrict",
            varargs_list = list(NULL), metadata_filter = NULL
        ),
        "umap;CTstrict;_;_"
    )

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap",
            clonecall = "CTstrict",
            varargs_list = list("seurat_clusters" = c(2, 5, 8)),
            metadata_filter = NULL
        ),
        "umap;CTstrict;seurat_clusters=c(2,5,8);_"
    )

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap",
            clonecall = "CTstrict",
            varargs_list = list("seurat_clusters" = c(8, 2, 5)),
            metadata_filter = NULL
        ),
        "umap;CTstrict;seurat_clusters=c(2,5,8);_"
    )

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap",
            clonecall = "CTstrict",
            varargs_list = list("orig.id" = c("P19B", "P19L", "P20B", "P20L")),
            metadata_filter = NULL
        ),
        "umap;CTstrict;orig.id=c('P19B','P19L','P20B','P20L');_"
    )

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap",
            clonecall = "CTstrict",
            varargs_list = list("orig.id" = c("P19L", "P20B", "P19B", "P20L")),
            metadata_filter = NULL
        ),
        "umap;CTstrict;orig.id=c('P19B','P19L','P20B','P20L');_"
    )

})

test_that("parse_to_object_id works for multiple varargs", {

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap",
            clonecall = "CTstrict",
            varargs_list = list(
                "seurat_clusters" = c(6, 2, 5),
                "orig.id" = c("P19L", "P20B", "P19B", "P20L")
            ),
            metadata_filter = NULL
        ),

        paste(
            "umap",
            "CTstrict",
            "orig.id=c('P19B','P19L','P20B','P20L'),seurat_clusters=c(2,5,6)",
            "_",
            sep = ";"
        )
    )

})