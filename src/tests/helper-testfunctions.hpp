#include <Rcpp.h>
#include <string>

bool approx_equal(double a, double b, double epsilon = 5e-5) {
    return std::abs(a - b) <= epsilon;
}

bool elements_are_equal(
    Rcpp::NumericVector vec1, Rcpp::NumericVector vec2, double epsilon = 5e-5
) {
    for (int i = 0; i < vec1.size(); i++) {
        if (!approx_equal(vec1[i], vec2[i], epsilon)) {
            return false;
        }
    }
    return true;
}