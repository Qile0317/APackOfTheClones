This folder is a work in progress directory of a c++ version of the `Main_algo.R` script to improve performance, and will be converted into an `Rcpp` API in the future.

for devs; currently the hierarchy of .h files are:
1. `data_structure.h`
2. `circle_packing.h`
And 1 is nested in 2. 2 is nested in testing.cpp. it should'nt be this way in the future.