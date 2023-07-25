context("C++ | utils") {
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

    test_that("elements_are_equal works") {
        Rcpp::NumericVector v1 = {1,2,3};
        Rcpp::NumericVector v2 = {1,2,3};
        expect_true(elements_are_equal(v1, v2));

        v2 = {1.1,2.1,3.1};
        expect_true(elements_are_equal(v1, v2, 0.15));
        expect_false(elements_are_equal(v1, v2));
    }

    test_that("has_repeats works") {
        Rcpp::NumericVector v1 = {1,2,3};
        Rcpp::NumericVector v2 = {4,1,5};
        expect_true(has_repeats(v1, v2));

        v1 = {1,2,3}; v2 = {4,6,5};
        expect_false(has_repeats(v1, v2));
    }

    test_that("has_common_strs works") {
        std::vector<std::string> v1 = {"a", "b", "c"};
        std::vector<std::string> v2 = {"d", "e", "c"};
        expect_true(has_common_strs(v1, v2));
    }
}
