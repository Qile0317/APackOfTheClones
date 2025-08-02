test_that("overlayLegend() and removeLegend() works", {
    data("combined_pbmc")

    expect_doppelganger(
        "overlayLegend plot",
        vizAPOTC(combined_pbmc, verbose = FALSE)
    )

    expect_doppelganger(
        "overlayLegend plot",
        vizAPOTC(combined_pbmc,  add_size_legend = FALSE, verbose = FALSE) %>%
            overlayLegend()
    )

    expect_doppelganger(
        "removeLegend plot",
        vizAPOTC(combined_pbmc, add_size_legend = FALSE, verbose = FALSE)
    )

    expect_doppelganger(
        "removeLegend plot",
        vizAPOTC(combined_pbmc, verbose = FALSE) %>% removeLegend()
    )

    expect_doppelganger(
        "overlayLegend plot",
        vizAPOTC(combined_pbmc, verbose = FALSE) %>%
            removeLegend() %>%
            overlayLegend()
    )

    expect_doppelganger(
        "removeLegend plot",
        vizAPOTC(combined_pbmc, verbose = FALSE) %>%
            removeLegend() %>%
            overlayLegend() %>%
            removeLegend()
    )
})
