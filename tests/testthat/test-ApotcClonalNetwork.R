quietly_test_that("getting shared clones works", {

    expect_contains(
        object = getSharedClones(
            get(data("combined_pbmc")),
            clonesize_range = c(1L, Inf)
        ),
        expected = getdata("ApotcClonalNetwork", "shared_clones")
    )

    expect_equal(
        object = getSharedClones(
            get(data("combined_pbmc")),
            clonesize_range = c(1L, Inf)
        ),
        expected = getdata("ApotcClonalNetwork", "shared_clones")
    )

})

test_that("adding shared clone links works", {
    expect_doppelganger(
        "default shared clone line link plot",
        vizAPOTC(
            get(data("combined_pbmc")),
            show_shared_clones = TRUE,
            verbose = FALSE
        )
    )
})