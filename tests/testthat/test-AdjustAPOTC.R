quietly_test_that("AdjustAPOTC works", {

  combined_pbmc <- RunAPOTC(get(data("combined_pbmc")), verbose = FALSE)
  combined_pbmc_apotc_obj <- getLastApotcData(combined_pbmc)

  expect_equal(
    combined_pbmc_apotc_obj,

    combined_pbmc %>%
      AdjustAPOTC(
        nudge_cluster = 1,
        nudge_vector = c(1L, 1L),
      ) %>%
      AdjustAPOTC(
        nudge_cluster = 1,
        nudge_vector = c(-1L, -1L),
      ) %>%
      getLastApotcData()
  )

  expect_equal(
    combined_pbmc_apotc_obj,

    combined_pbmc %>%
      AdjustAPOTC(
        nudge_cluster = 1:17,
        nudge_vector = c(1L, 1L),
      ) %>%
      AdjustAPOTC(
        nudge_cluster = 1:17,
        nudge_vector = c(-1L, -1L),
      ) %>%
      getLastApotcData()
  )

  expect_equal(
    combined_pbmc_apotc_obj,

    combined_pbmc %>%
      AdjustAPOTC(
        nudge_cluster = 1:17,
        nudge_vector = init_list(17, c(831L, 317L)),
      ) %>%
      AdjustAPOTC(
        nudge_cluster = 1:17,
        nudge_vector = init_list(17, c(-831L, -317L)),
      ) %>%
      getLastApotcData()
  )

})

quietly_test_that("AdjustAPOTC doesn't affect empty clusters", {

  combined_pbmc <- RunAPOTC(
    get(data("combined_pbmc")), override = TRUE, verbose = FALSE
  )
  combined_pbmc_apotc_obj <- getLastApotcData(combined_pbmc)

  # move and nudge all empty clusters by random amounts
  # only the last cluster in the testdata is empty

  empty_cluster_inds <- get_clusterlists(combined_pbmc_apotc_obj) %>%
    sapply(is_empty) %>%
    which()

  set.seed(831)

  # test equality
  expect_equal(
    combined_pbmc_apotc_obj,

    combined_pbmc %>%
      AdjustAPOTC(
        relocate_cluster = empty_cluster_inds,
        relocation_coord = sample(-100:100, 2, replace = TRUE)
      ) %>%
      AdjustAPOTC(
        nudge_cluster = empty_cluster_inds,
        nudge_vector = sample(-100:100, 2, replace = TRUE)
      ) %>%
      getLastApotcData()
  )

})
