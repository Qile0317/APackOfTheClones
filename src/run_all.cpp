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

// function scripts
#include "get_clone_sizes.h"
#include "circle_layout.h"
#include "repulsion.h"

// test scripts for testthat
#include "test-circle_layout.h"