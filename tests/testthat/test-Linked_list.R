test_that("Node initialization works", {
  test_node <- node$new(
    val = list(1, 1, 2, 3, "color" = "#5318008", 3)
  )

  expect_true(is.R6(test_node))
  expect_equal(test_node$val[[2]], 1)
  expect_equal(test_node$val[[4]], 3)
})
