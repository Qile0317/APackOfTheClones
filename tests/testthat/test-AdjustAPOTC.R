quietly_test_that("AdjustAPOTC nudging works", {

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
        nudge_label = 1:17,
        label_nudge_vector = init_list(17, c(831L, 317L))
      ) %>%
      AdjustAPOTC(
        nudge_cluster = 1:17,
        nudge_vector = init_list(17, c(-831L, -317L)),
        nudge_label = 1:17,
        label_nudge_vector = init_list(17, c(-831L, -317L))
      ) %>%
      getLastApotcData()
  )

})

# not comprehensive at all
quietly_test_that("AdjustAPOTC's other args work", {

  combined_pbmc <- RunAPOTC(get(data("combined_pbmc")), override = TRUE)

  expect_equal(
    expected = combined_pbmc,
    object = AdjustAPOTC(
      combined_pbmc,
      new_rad_scale_factor = 0.95,
      recolor_cluster = 1:17,
      new_color = c(
        "#F8766D", "#E7851E", "#D09400", "#B2A100", "#89AC00", "#45B500",
        "#00BC51", "#00C087", "#00C0B2", "#00BCD6", "#00B3F2", "#29A3FF",
        "#9C8DFF", "#D277FF", "#F166E8", "#FF61C7", "#FF689E"
      ),
      rename_label = 1:17,
      new_label = c(
        "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
        "C11", "C12", "C13", "C14", "C15", "C16", "C17"
      )
    )
  )

  combined_pbmc %>% expect_equal(expected = .,
    object = AdjustAPOTC(
      .,
      new_clone_scale_factor = get_clone_scale_factor(getLastApotcData(.))
    )
  )

})

quietly_test_that("AdjustAPOTC doesn't affect empty clusters", {

  data("combined_pbmc")
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
