quietly_test_that("getting shared clones works", {

    data("combined_pbmc")
    expected_all_shared <- getdata("ApotcClonalNetwork", "shared_clones")

    expect_equal(
        object = getSharedClones(combined_pbmc),
        expected = expected_all_shared
    )

    # test top

    expect_equal(
        object = getSharedClones(combined_pbmc, top = 1L),
        expected = expected_all_shared[2]
    ) # TODO verify

    expect_equal(
        object = getSharedClones(combined_pbmc, top = 1 / 9),
        expected = expected_all_shared[2]
    )

    expect_equal(
        object = getSharedClones(combined_pbmc, top = 10L),
        expected = expected_all_shared
    )

    # check everything = Inf does nothing

    expect_equal(
        object = getSharedClones(
            combined_pbmc,
            top = Inf,
            top_per_cl = Inf,
            intop = Inf,
            intop_per_cl = Inf
        ),
        expected = expected_all_shared
    )

    # check empty shared clones

    # TODO everything else needs verification

})

# TODO following two tests unfinished

test_that("filter_top_clones works", {
    clone_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")
    expect_equal(
        filter_top_clones(clone_sizes, Inf),
        clone_sizes
    )
})

test_that("filter_top_by_cluster works", {
    clone_sizes <- getdata("get_clone_sizes", "raw_strict_clone_sizes")
    expect_equal(
        object = filter_top_by_cluster(clone_sizes, Inf),
        expected = sort_each_table(clone_sizes, desc = TRUE)
    )
})

test_that("adding shared clone links works", {
    data("combined_pbmc")
    expect_doppelganger(
        "default shared clone line link plot",
        vizAPOTC(
            combined_pbmc,
            show_shared = getSharedClones(combined_pbmc),
            seurat_clusters = c(3:5, 9, 13),
            add_size_legend = FALSE,
            retain_axis_scales = TRUE,
            verbose = FALSE
        )
    )
})
