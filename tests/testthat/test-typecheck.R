test_that("basic typechecking works", {

  a <- 10
  y <- list(a = 10)
  z <- rep(8, 31)

  expect_null(typecheck(a, is_an_integer, is_integer, is_a_numeric))

  expect_error(
    object = typecheck(a, is_a_character),
    "`a` must be a character of length 1"
  )

  expect_error(
    object = typecheck(y, is_an_integer),
    "`y` must be an integer of length 1"
  )

  expect_null(typecheck(y$a, is_integer))

  expect_error(
    object = typecheck(y$a, is_a_character),
    "`a` must be a character of length 1"
  )

  expect_error(
    object = typecheck(z, is_a_logical, is_character),
    "`z` must be a character vector, or a logical of length 1"
  )

})
