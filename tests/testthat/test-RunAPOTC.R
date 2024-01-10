# TODO RunAPOTC tests

quietly_test_that("exception handling in RunAPOTC works", {
    data("combined_pbmc")

    # test override
    expect_error(
        RunAPOTC(RunAPOTC(combined_pbmc)),
        paste(
            "An APackOfTheClones run with the the parameters umap;CTstrict;_;_",
            "appears to already have been ran. If this is intended, set the",
            "`override` argument to `TRUE` and re-run."
        ),
        fixed = TRUE
    )

    # test wrong arguments
    # - due to how R works internally with ...,
    #   prefixes of the first 3 args wont work!
    #   E.g. `reduction`/`clone` actually works as well!
    expect_error(
        RunAPOTC(combined_pbmc, ID = "r1"),
        c(paste(
            "should `ID` be a function argument? ",
            "If so, did you mean `run_id`? ",
            "Otherwise, should `ID` be a subsetting argument? ",
            "If so, did you mean `CTnt`?", sep = ""
        )),
        fixed = TRUE
    )

    # TODO more tests

})
