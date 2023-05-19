// this is the file to execute all scripts and test scripts to handle both 
// exported and internal cpp functions

// header deps
#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

// convenience functions and definitions
using Llui = long long unsigned int;

inline double sqr(double n) {
  return n*n;
}

bool approx_equal(double a, double b, double epsilon = 5e-5) {
  return std::abs(a - b) <= epsilon;
}

// function scripts
#include "get_clone_sizes.h"
#include "circle_layout.h"
#include "repulsion.h"

// test scripts for testthat
#include "tests/test-circle_layout.h"