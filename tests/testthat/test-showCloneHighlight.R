# test_that("showCloneHighlight works", {

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

#     # TODO test S4 also present in scRepertoire

#     # TODO test missing clonotypes
# })
