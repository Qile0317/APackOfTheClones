# test_that("vdiffr test of showCloneHighlight works", {

#     skip_on_ci() # there is no visually distinguishable difference in the artifact

#     test_cl_5_9_highlighted_ggplot <- get(data("combined_pbmc")) %>%
#         vizAPOTC(
#             clonecall = "aa",
#             seurat_clusters = c(5, 9),
#             retain_axis_scales = TRUE,
#             add_size_legend = FALSE,
#             verbose = FALSE
#         ) %>%
#         showCloneHighlight(c(
#             "CASLSGSARQLTF_CASSPTVAGEQFF", "CAELNQAGTALIF_CASSQAPFSTSGELFF"
#         ))

#     expect_doppelganger(
#         "showCloneHighlight_only_clusters_5_9",
#         test_cl_5_9_highlighted_ggplot
#     )

# })

test_that("showCloneHighlight works", {

    data("combined_pbmc")

    test_cl_5_9_highlighted_data <- combined_pbmc %>%
        vizAPOTC(
            clonecall = "aa",
            seurat_clusters = c(5, 9),
            retain_axis_scales = TRUE,
            add_size_legend = FALSE,
            verbose = FALSE
        ) %>%
        showCloneHighlight(c(
            "CASLSGSARQLTF_CASSPTVAGEQFF", "CAELNQAGTALIF_CASSQAPFSTSGELFF"
        )) %>%
        get_ggplot_data()

    expect_equal(
        object = test_cl_5_9_highlighted_data,
        expected = getdata(
            "showCloneHighlight", "cl_5_9_highlighted_data"
        )
    )

    # TODO test S4 also present in scRepertoire

    # TODO test missing clonotypes

    expect_warning(
        vizAPOTC(
            combined_pbmc, clonecall = "aa", verbose = FALSE
        ) %>% showCloneHighlight(
            c("CASLSGSARQLTF_CASSPTVAGEQFF", "CAELNQAGTALIF_CASSQAPFSTSGELFF"),
            color_each = FALSE, default_color = NULL, scale_bg = 1
        )
    )

})