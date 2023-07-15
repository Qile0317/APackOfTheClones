// this is the file to execute all scripts and test scripts to handle both 
// exported and internal cpp functions

// dependency headers
#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <testthat.h>

// this can be removed when the C++ code is fixed
// [[Rcpp::export]]
double estimate_rad(
    std::vector<double> x_vals,
    std::vector<double> rad_vals,
    double centroid_x
) {
  double max_x = 0;
  int max_ind = 0, n = x_vals.size();
  for (int i = 0; i < n; i++) {
    if (x_vals[i] > max_x) {
      max_x = x_vals[i];
      max_ind += i - max_ind;
    }
  }
  double result_num = max_x + rad_vals[max_ind] - centroid_x;
  return result_num;
}

// function scripts
#include "utils.h"
#include "get_clone_sizes.h"
#include "cpp_circle_layout.h"
#include "repulsion.h"

// test scripts for testthat
#include "tests/test-utils.hpp"
#include "tests/test-cpp_circle_layout.hpp"
#include "tests/test-repulsion.hpp"