test_that("parse_to_metadata_filter_str works for no varargs", {

    expect_identical(parse_to_metadata_filter_str(NULL, list()), "")
    expect_identical(parse_to_metadata_filter_str("", list()), "")
    expect_identical(parse_to_metadata_filter_str("   ", list()), "")

    expect_identical(
        parse_to_metadata_filter_str(
            "mito.genes > 1.5 | length(CTaa) - 1 < 20", list()
        ),
        "mito.genes>1.5|length(CTaa)-1<20"
    )
})

test_that("parse_to_metadata_filter_str works for base cases", {

    expect_identical(
        parse_to_metadata_filter_str(
            NULL, list("seurat_clusters" = 1)
        ),
        "seurat_clusters==1"
    )

    expect_identical(
        parse_to_metadata_filter_str(
            NULL, list("seurat_clusters" = c(2, 5, 8))
        ),
        "seurat_clusters==2|seurat_clusters==5|seurat_clusters==8"
    )

    expect_identical(
        parse_to_metadata_filter_str(
            NULL, list("seurat_clusters" = c(8, 2, 2, 5, 2))
        ),
        "seurat_clusters==2|seurat_clusters==5|seurat_clusters==8"
    )

    expect_identical(
        parse_to_metadata_filter_str(
            NULL, list("orig.id" = "P19B")
        ),
        "orig.id=='P19B'"
    )

    expect_identical(
        parse_to_metadata_filter_str(
            NULL, list("orig.id" = c("P19B", "P19L"))
        ),
        "orig.id=='P19B'|orig.id=='P19L'"
    )

    expect_identical(
        parse_to_metadata_filter_str(
            NULL, list("orig.id" = c("P19L", "P19B", "P19B"))
        ),
        "orig.id=='P19B'|orig.id=='P19L'"
    )

})

test_that("parse_to_metadata_filter_str works with all args", {

    expect_identical(
        parse_to_metadata_filter_str(
            metadata_filter = "mito.genes > 1.5 | length(CTaa) - 1 < 20",
            varargs_list = list(
                "seurat_clusters" = c(2, 5, 8),
                "orig.id" = c("P19B", "P19L")
            )
        ),
        paste(
            "((orig.id=='P19B'|orig.id=='P19L')",
            "(seurat_clusters==2|seurat_clusters==5|seurat_clusters==8))",
            "(mito.genes>1.5|length(CTaa)-1<20)",
            sep = "&"
        )
    )

    expect_identical(
        parse_to_metadata_filter_str(
            metadata_filter = "mito.genes > 1.5 | length(CTaa) - 1 < 20",
            varargs_list = list(
                "seurat_clusters" = c(8, 2, 2, 5, 2),
                "orig.id" = c("P19L", "P19B", "P19B")
            )
        ),
        paste(
            "((orig.id=='P19B'|orig.id=='P19L')",
            "(seurat_clusters==2|seurat_clusters==5|seurat_clusters==8))",
            "(mito.genes>1.5|length(CTaa)-1<20)",
            sep = "&"
        )
    )

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

    expect_identical(
        parse_to_object_id(
            reduction_base = "umap",
            clonecall = "CTstrict",
            varargs_list = list(
                "seurat_clusters" = c(6, 2, 0, 2, 5),
                "orig.id" = c("P20B", "P19L", "P20B", "P19B", "P20L", "P20B")
            ),
            metadata_filter = NULL
        ),

        paste(
            "umap",
            "CTstrict",
            "orig.id=c('P19B','P19L','P20B','P20L'),seurat_clusters=c(0,2,5,6)",
            "_",
            sep = ";"
        )
    )

})

test_that("parse_to_object_id works with all args changed", {

    expect_identical(
        parse_to_object_id(
            reduction_base = "tsne",
            clonecall = "CTnt",
            varargs_list = list(
                "seurat_clusters" = c(6, 2, 0, 2, 5),
                "orig.id" = c("P20B", "P19L", "P20B", "P19B", "P20L", "P20B")
            ),
            metadata_filter = "mito.genes > 1.5 | length(CTaa) - 1 < 20"
        ),

        paste(
            "tsne",
            "CTnt",
            "orig.id=c('P19B','P19L','P20B','P20L'),seurat_clusters=c(0,2,5,6)",
            "mito.genes>1.5|length(CTaa)-1<20",
            sep = ";"
        )
    )

})
