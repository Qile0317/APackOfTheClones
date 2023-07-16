// interesting property: calling .size() on a NumericVector with NULLs
// will ignore the null values

//function to transform the raw clone sizes for circle_layout
// [[Rcpp::export]]
Rcpp::List get_transformed_clone_sizes(
  Rcpp::List sizelist, double clone_scale_factor, int num_clusters
) {
  Rcpp::List output_sizes (num_clusters);
  for (int i = 0; i < num_clusters; i++) {
    int n = sizelist[i].size();
    Rvec curr_sizes = sizelist[i];
    if (curr_sizes.isNotNull()) {
      for (int j = 0; j < n; j++) {
        curr_sizes[j] = sqrt(curr_sizes[j]) * clone_scale_factor;
      }
      output_sizes[i] = curr_sizes;
    } else {
      output_sizes[i] = Rcpp::List::create(); // if nothing, creates empty list
    }
  }
  return output_sizes;
}

// tested in the R testthat folder