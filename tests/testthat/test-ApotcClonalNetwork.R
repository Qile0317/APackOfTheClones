quietly_test_that("getting shared clones works", {

    expect_equal(
        object = getSharedClones(get(data("combined_pbmc"))),
        expected = getdata("ApotcClonalNetwork", "shared_clones")
    )

})

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
    expect_doppelganger(
        "default shared clone line link plot",
        vizAPOTC(
            get(data("combined_pbmc")),
            show_shared_clones = TRUE,
            seurat_clusters = c(3:5, 9, 13),
            add_size_legend = FALSE,
            retain_axis_scales = TRUE,
            verbose = FALSE
        )
    )
})
