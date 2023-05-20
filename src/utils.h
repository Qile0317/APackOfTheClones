using Llui = long long unsigned int;
using Rcpp::_; // for overlap_check
using Rvec = Rcpp::Nullable<Rcpp::NumericVector>; // for repulsion

#define all(x) (x).begin(), (x).end()

inline double sqr(double n) {
  return n*n;
}

bool approx_equal(double a, double b, double epsilon = 5e-5) {
  return std::abs(a - b) <= epsilon;
}

// copy of the R progress_bar. not unit tested but works perfectly
void progress_bar(int x, int max) {
  double percent = 100.0 * (double(x) / double(max));
  int filledWidth = static_cast<int>(percent * 0.5);
  int remainingWidth = 50 - filledWidth;
  
  std::string progressBar = "\r[";
  progressBar += std::string(filledWidth, '=');
  progressBar += std::string(remainingWidth, ' ') + "] ";
  progressBar += std::to_string(int(floor(percent))) + "%";
  
  Rcpp::Rcout << progressBar;
}