context("Cpp utils") {
  test_that("sqr works") {
    expect_true(sqr(-5) == 25);
    expect_true(sqr(12) == 144);
  }
  
  test_that("appox_equal works") {
    expect_true(approx_equal(1, 1));
    expect_true(approx_equal(1.0001, 1, 1e-4));
    expect_false(approx_equal(2, 1));
    expect_false(approx_equal(1.0001, 1, 1e-5));
  }
}
