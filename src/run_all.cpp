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

// function scripts
#include "utils.h"
#include "get_clone_sizes.h"
#include "cpp_circle_layout.h"
#include "repulsion.h"

// test scripts for testthat
#include "tests/test-utils.hpp"
#include "tests/test-cpp_circle_layout.hpp"
#include "tests/test-repulsion.hpp"