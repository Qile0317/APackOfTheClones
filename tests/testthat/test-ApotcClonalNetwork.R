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

    # expect_identical(
    #     getSharedClones(
    #         get(data("combined_pbmc")),
    #         clonesize_range = c(100L, Inf)
    #     ),
    #     list()
    # )
})
