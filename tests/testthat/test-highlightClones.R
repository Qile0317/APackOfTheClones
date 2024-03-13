test_that("highlightClones works", {
    expect_doppelganger(
        "highlightClones_only_clusters_5_9",
        get(data("combined_pbmc")) %>%
            vizAPOTC(
                clonecall = "aa",
                seurat_clusters = c(5, 9),
                retain_axis_scales = TRUE,
                add_size_legend = FALSE,
                verbose = FALSE
            ) %>%
            APackOfTheClones::highlightClones(c(
                "CASLSGSARQLTF_CASSPTVAGEQFF", "CAELNQAGTALIF_CASSQAPFSTSGELFF"
            ))
    )

    # TODO test S4 also present in scRepertoire

    # TODO test missing clonotypes
})