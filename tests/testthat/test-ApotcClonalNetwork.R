test_that("getting shared clones works", {

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

    # TODO test getting top 4 - the first should be the same (size 11) while 2,3,4, (size 3) should setequal

    expect_equal(
        object = getSharedClones(combined_pbmc, top = 10L),
        expected = expected_all_shared
    )

    # test intop

    expect_identical(
        object = getSharedClones(combined_pbmc, intop = 1),
        expected = expected_all_shared[2]
    )

    expect_identical(
        object = getSharedClones(combined_pbmc, intop = 2),
        expected = expected_all_shared[c(2, 5, 7)]
    )

    expect_mapequal(
        object = getSharedClones(combined_pbmc, intop = 3),
        expected = expected_all_shared
    )

    expect_mapequal(
        object = getSharedClones(combined_pbmc, intop = Inf),
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

    # check empty shared clones if seurat object is filtered

    expect_identical(
        object = getSharedClones(combined_pbmc, seurat_clusters = 15:17),
        expected = list()
    )

    expect_identical(
        object = getSharedClones(combined_pbmc, seurat_clusters = 12),
        expected = list()
    )

    expect_identical(
        object = getSharedClones(combined_pbmc, publicity = c(3, Inf)),
        expected = expected_all_shared[5]
    )

    # TODO check the intersection
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
            # also a test of these two args
            retain_axis_scales = TRUE,
            show_labels = TRUE,
            verbose = FALSE
        )
    )
})
