/*
#include <Rcpp.h>
#include <vector>
#include <cmath>
*/
 
using Rvec = Rcpp::Nullable<Rcpp::NumericVector>;
using Llui = long long unsigned int;

// interesting property: calling .size() on a NumericVector with NULLs
// will ignore the null values

//function to transform the raw clone sizes for circle_layout
// [[Rcpp::export]]
Rcpp::List get_transformed_clone_sizes(
  Rcpp::List sizelist, double clone_scale_factor, int num_clusters
) {
  Rcpp::List output_sizes = Rcpp::List::create();
  for (int i = 0; i < num_clusters; i++) {
    Rvec curr_sizes = sizelist[i];
    if (curr_sizes.isNotNull()) {
      std::vector<double> curr_sizes (curr_sizes); 
      for (Llui j = 0; j < curr_sizes.size(); j++) {
        curr_sizes[j] = sqrt(curr_sizes[j]) * clone_scale_factor;
      }
      output_sizes[i] = curr_sizes;
    } else {
      output_sizes[i] = Rcpp::List::create(); // if nothing, creates empty list
    }
  }
  return output_sizes;
}