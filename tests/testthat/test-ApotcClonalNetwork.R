quietly_test_that("getting shared clones works", {
    expect_identical(
        getSharedClones(
            get(data("combined_pbmc")),
            clonesize_range = c(1L, Inf)
        ),
        getdata("ApotcClonalNetwork", "shared_clones")
    )

    # expect_identical(
    #     getSharedClones(
    #         get(data("combined_pbmc")),
    #         clonesize_range = c(100L, Inf)
    #     ),
    #     list()
    # )
})
