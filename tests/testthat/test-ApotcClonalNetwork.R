test_that("getting shared clones works", {
    expect_identical(
        getSharedClones(get(data("combined_pbmc"))),
        getdata("ApotcClonalNetwork", "shared_clones")
    )
})
