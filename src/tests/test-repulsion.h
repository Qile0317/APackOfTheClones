context("Cpp repulsion functions") {
  test_that("neg_dir_vec works") {
    expect_true(1==1); // unfinished
  }
  
  test_that("polar_form works") {
    std::vector<double> trial = polar_form({1, -2});
    expect_true(approx_equal(trial[0], sqrt(5)));
    expect_true(approx_equal(trial[1], -1.10714872));
  }
  
  test_that("neg_polar_dist_vec works") {
    expect_true(1==1);  // unfinished
  }
  
  test_that("component_form works") {
    std::vector<double> trial = component_form({sqrt(5), -1.10714872});
    expect_true(approx_equal(trial[0], 1));
    expect_true(approx_equal(trial[1], -2));
  }
  
  // all Rcpp exported functions are tested in R
}
