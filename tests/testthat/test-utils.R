test_that("getReductionCentroids works", {
    expect_equal(
        getReductionCentroids(get(data("combined_pbmc")), "umap"),
        getdata("combined_pbmc", "all_cluster_centroids"),
        tolerance = 1e-6
    )
})

test_that("progress_bar works", {
    expect_identical(
        capture.output(progress_bar(1, 1)),
        "\r[==================================================] 100%"
    )
    expect_identical(
        capture.output(progress_bar(69, 420)),
        "\r[========                                          ] 16%"
    )
})

test_that("start_progress_bar works", {
    expect_identical(
        capture.output(start_progress_bar(TRUE)),
        "\r[                                                  ] 0%"
    )
})

test_that("start_progress_bar works", {
    expect_identical(
        capture.output(end_progress_bar(TRUE)),
        "\r[==================================================] 100%"
    )
})

# test table utils

test_that("convert_named_numeric_to_table works", {

    expect_identical(
        object = convert_named_numeric_to_table(c("b" = 10, "c" = 12)),
        expected = table(rep(c("b", "c"), times = c(10, 12)))
    )

    expect_identical(
        object = convert_named_numeric_to_table(c("b" = 10)),
        expected = table(rep("b", times = 10))
    )

})

test_that("union_list_of_tables works", {

    expect_identical(
        union_list_of_tables(list(table(letters[1:5]), table(letters[1:5]))),
        expected = strip_to_numeric(table(rep(letters[1:5], 2)))
    )

    expect_identical(
        union_list_of_tables(list(table(letters))),
        expected = init_numeric(letters, 1)
    )

    expect_identical(
        union_list_of_tables(init_list(5, table(letters))),
        expected = init_numeric(letters, 5)
    )

})

# test readability

test_that("isnt_empty works", {
    expect_true(isnt_empty(list(c(1, 2, 3))))
    expect_false(isnt_empty(list()))
})

test_that("isnt_na works", {
    expect_true(isnt_na(list(c(1, 2, 3))))
    expect_false(isnt_na(list(NA, NA)))
})

test_that("isnt_empty_nor_na works", {
    expect_true(isnt_empty_nor_na(list(c(1, 2, 3))))
    expect_false(isnt_empty_nor_na(list(NA, NA)))
    expect_false(isnt_empty_nor_na(list()))
})

test_that("is_int works", {
    expect_true(is_int(1L))
    expect_false(is_int(1.001))
})

test_that("should_estimate works", {
    expect_true(should_estimate("auto"))
    expect_false(should_estimate(10290987))
})

test_that("get_xr and get_yr works", {
    plt <- ggplot2::ggplot(data.frame(x = c(1, 2), y = c(3, 4))) +
        ggplot2::geom_point(ggplot2::aes(x, y))

    expect_equal(get_xr(plt), c(1, 2))
    expect_equal(get_yr(plt), c(3, 4))
    # TODO more tests
})

test_that("prepend_indefinite_article works", {
    expect_identical(prepend_indefinite_article("A"), "an A")
    expect_identical(prepend_indefinite_article("B"), "a B")
    expect_identical(prepend_indefinite_article("NULL"), "NULL")
})

test_that("subsetSeuratMetaData works", {
    # TODO more tests

    expect_error(
        subsetSeuratMetaData(get(data("combined_pbmc")), "seurat_clusters==18"),
        "please check `extra_filter`, no rows in the seurat metadata match the filter condition"
    )

})

quietly_test_that("attempt_correction works", {
    data("combined_pbmc")

    combined_pbmc@reductions[["pca"]] <- 0
    expect_identical("umap", attempt_correction(combined_pbmc, "Umap"))
    expect_identical("pca", attempt_correction(combined_pbmc, "PCA"))

    combined_pbmc@reductions[["tsne"]] <- 0
    expect_identical("tsne", attempt_correction(combined_pbmc, "t-SNE"))

    combined_pbmc@reductions[["fooBarReduction"]] <- 0
    combined_pbmc@reductions[["Trex.umap"]] <- 0
    expect_identical(
        "fooBarReduction", attempt_correction(combined_pbmc, "foobarreduction")
    )
    expect_identical(
        "Trex.umap", attempt_correction(combined_pbmc, " trex. umap")
    )

    combined_pbmc@reductions[["Foo 'bar' baz"]] <- 0
    combined_pbmc@reductions[["Foo ' b  a r ' baz"]] <- 0
    expect_identical(
        "Foo 'bar' baz", attempt_correction(combined_pbmc, "foo'bar'baz")
    )
    expect_identical(
        "Foo ' b  a r ' baz",
        attempt_correction(combined_pbmc, "foo' b  a r 'baz")
    )

    # TODO more tests
})

test_that("closest_word works", {
    expect_identical(closest_word(" umsp", c("umap", "tsne", "pca")), "umap")
    expect_identical(closest_word("t snq", c("umap", "tsne", "pca")), "tsne")
    expect_identical(closest_word("aca ", c("umap", "tsne", "pca")), "pca")
})

test_that("construct_prefix_vector works", {
    test_samples <- c("NP1", "NP2", "NP1", "NP2")
    test_ids <- c("CTRL", "CTRL", "STIM", "STIM")

    test_results <- construct_prefix_vector(list(test_samples, test_ids))
    expect_identical(
        test_results, c("NP1_CTRL", "NP2_CTRL", "NP1_STIM", "NP2_STIM")
    )

    test_age <- c("Y", "O", "Y", "O")

    test_results <- construct_prefix_vector(
        list(test_age, test_samples, test_ids)
    )

    expect_identical(
        test_results, c("Y_NP1_CTRL", "O_NP2_CTRL", "Y_NP1_STIM", "O_NP2_STIM")
    )
})

# test math utils

test_that("get_unique_pairs_up_to works", {

    expect_identical(get_unique_pairs_up_to(1), list())

    expect_identical(get_unique_pairs_up_to(2), list(1:2))

    expect_identical(
        get_unique_pairs_up_to(5),
        list(
            1:2, c(1L, 3L), c(1L, 4L), c(1L, 5L), 2:3, c(2L, 4L), c(2L, 5L),
            3:4, c(3L, 5L), 4:5
        )
    )
})

# test spelling utils

test_that("strip_unquoted_spaces works", {
    expect_identical(strip_unquoted_spaces("foo"), "foo")
    expect_identical(strip_unquoted_spaces(c("foo", "bar")), c("foo", "bar"))

    expect_identical(strip_unquoted_spaces(" foo "), "foo")
    expect_identical(
        strip_unquoted_spaces(c(" foo ", " bar ")), c("foo", "bar")
    )

    expect_identical(strip_unquoted_spaces(" f oo "), "foo")
    expect_identical(
        strip_unquoted_spaces(c(" f oo ", " b ar ")), c("foo", "bar")
    )

    expect_identical(strip_unquoted_spaces(" f f ' o o ' "), "ff' o o '")
    expect_identical(
        strip_unquoted_spaces(c(" f f ' o o ' ", " bb b ' a r ' rr r ")),
        c("ff' o o '", "bbb' a r 'rrr")
    )
})

test_that("get_ident_levels works", {
    data("combined_pbmc")

    expect_identical(
        object = get_ident_levels(combined_pbmc, custom_ident = NULL),
        expected = get_idents(ApotcData(
            seurat_obj = combined_pbmc,
            alt_ident = NULL,
            metadata_filter_condition = "",
            clonecall = "CTstrict",
            reduction_base = "umap",
            clone_scale_factor = 0.300051,
            rad_scale_factor = 0.95
        ))
    )

    expect_identical(
        object = get_ident_levels(combined_pbmc, custom_ident = "orig.ident"),
        expected = get_idents(ApotcData(
            seurat_obj = combined_pbmc,
            alt_ident = "orig.ident",
            metadata_filter_condition = "",
            clonecall = "CTstrict",
            reduction_base = "umap",
            clone_scale_factor = 0.300051,
            rad_scale_factor = 0.95
        ))
    )

})

# TODO more testcases