expect_apotc_data_clusterlist_equals_expected <- function(
    apotc_obj, expected_clusterlist, tolerance = 1e-6
) {
    testthat::expect_true(
        does_apotc_data_clusterlist_equal_expected(
            apotc_obj, expected_clusterlist, tolerance, verbose = TRUE
        )
    )
}
