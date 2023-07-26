data("mini_seurat_obj", "mini_clonotype_data")

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

# print_completion_time cant really be tested but is extremely simple

test_that("isnt_empty works", {
    expect_true(isnt_empty(list(c(1,2,3))))
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
})

test_that("attempt_correction works", {
    expect_identical("umap", attempt_correction("Umap"))
    expect_identical("tsne", attempt_correction("t-SNE"))
    expect_identical("pca", attempt_correction("PCA"))
})

test_that("closest_word works", {
    expect_identical(closest_word(" umsp"), "umap")
    expect_identical(closest_word("t snq"), "tsne")
    expect_identical(closest_word("aca "), "pca")
})

test_that("metadata_name_warnstring works", {
    expect_null(metadata_name_warnstring(mini_seurat_obj, mini_clonotype_data))

    fake_tcr_df <- data.frame(c(1,2),c(3,4))
    names(fake_tcr_df) <- c("NA", NA)

    expect_identical(
        metadata_name_warnstring(mini_seurat_obj, fake_tcr_df),
        "tcr_dataframe has NAs in names, please fix"
    )

    names(fake_tcr_df) <- c("NA", "orig.ident")

    expect_identical(
        metadata_name_warnstring(mini_seurat_obj, fake_tcr_df),
        "tcr_dataframe has repeated names with the seurat_object metadata"
    )
})
