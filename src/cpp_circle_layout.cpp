#include <Rcpp.h>
#include <vector>
#include "CirclePacker.h"

// [[Rcpp::export]]
Rcpp::List cpp_circle_layout(
    std::vector<double> input_rad_vec,
    Rcpp::NumericVector centroid,
    double rad_decrease = 0,
    bool try_place = false,
    bool verbose = true
) {
    return CirclePacker::pack(input_rad_vec, centroid, rad_decrease, try_place, verbose);
}
