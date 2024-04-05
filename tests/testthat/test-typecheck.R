test_that("basic typechecking works", {

  x <- 10
  y <- list(x = 10)
  z <- rep(8, 31)

  expect_null(typecheck(x, is_an_integer, is_integer, is_a_numeric))

  expect_error(
    object = typecheck(x, is_a_character),
    "`x` must be a character of length 1"
  )

  expect_error(
    object = typecheck(y, is_an_integer),
    "`y` must be an integer of length 1"
  )

  expect_null(typecheck(y$x, is_integer))

  expect_error(
    object = typecheck(y$x, is_a_character),
    "`x` must be a character of length 1"
  )

  expect_error(
    object = typecheck(z, is_a_logical, is_character),
    "`z` must be a character vector, or a logical of length 1"
  )

})
