runTestShowCloneHighlight <- function() {
    get(data("combined_pbmc")) %>%
        vizAPOTC(
            clonecall = "aa",
            seurat_clusters = c(5, 9),
            retain_axis_scales = TRUE,
            add_size_legend = FALSE,
            verbose = FALSE
        ) %>%
        showCloneHighlight(c(
            "CASLSGSARQLTF_CASSPTVAGEQFF", "CAELNQAGTALIF_CASSQAPFSTSGELFF"
        ))
}

test_that("vdiffr test of showCloneHighlight works", {

    #skip_on_ci() # there is no visually distinguishable difference in the artifact

    expect_doppelganger(
        "showCloneHighlight_only_clusters_5_9",
        runTestShowCloneHighlight()
    )

})

test_that("showCloneHighlight works", {

    test_cl_5_9_highlighted_data <-
        get_ggplot_data(runTestShowCloneHighlight())

    # check the cluster, clonotype (and its color), and corresponding radii
    # are setequal since order doesn't matter
    expect_equal(
        object = test_cl_5_9_highlighted_data %>%
            dplyr::select(label, color, clonotype, r) %>%
            dplyr::arrange(clonotype),
        expected = getdata(
            "showCloneHighlight", "cl_5_9_highlighted_data"
        ) %>%
            dplyr::select(label, color, clonotype, r) %>%
            dplyr::arrange(clonotype)
    )

    # check that for each cluster and radii sets, there are setequal (x, y)
    # pairs since it would be geometrically equal
    expect_equal(
        object = test_cl_5_9_highlighted_data %>%
            dplyr::select(label, x, y, r) %>%
            dplyr::arrange(x, y),
        expected =
            getdata(
                "showCloneHighlight", "cl_5_9_highlighted_data"
            ) %>%
            dplyr::select(label, x, y, r) %>%
            dplyr::arrange(x, y),
        tolerance = 1e-5
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
