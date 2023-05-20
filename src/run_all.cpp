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
#include "circle_layout.h"
#include "repulsion.h"

// test scripts for testthat
#include "tests/test-utils.h"
#include "tests/test-circle_layout.h"
#include "tests/test-repulsion.h"