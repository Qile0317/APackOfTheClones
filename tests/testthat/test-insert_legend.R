test_that("overlayLegend() and removeLegend() works", {
    data("combined_pbmc")

    expect_equal(
        vizAPOTC(combined_pbmc, verbose = FALSE),
        vizAPOTC(combined_pbmc, verbose = FALSE) %>% overlayLegend()
    )

    expect_equal(
        vizAPOTC(combined_pbmc, add_size_legend = FALSE, verbose = FALSE),
        vizAPOTC(combined_pbmc, verbose = FALSE) %>% removeLegend()
    )

    expect_equal(
        vizAPOTC(combined_pbmc, verbose = FALSE),
        vizAPOTC(combined_pbmc, verbose = FALSE) %>%
            removeLegend() %>%
            overlayLegend()
    )
})
