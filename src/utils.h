using Llui = long long unsigned int;
using Rcpp::_; // for overlap_check

#define all(x) (x).begin(), (x).end()

inline double sqr(double n) {
  return n*n;
}

bool approx_equal(double a, double b, double epsilon = 5e-5) {
  return std::abs(a - b) <= epsilon;
}

void progress_bar(int x, int max) {
  double percent = 100.0 * (double(x) / double(max));
  int filledWidth = static_cast<int>(percent * 0.5);
  int remainingWidth = 50 - filledWidth;
  
  std::string progressBar;
  progressBar += "\r[";
  progressBar += std::string(filledWidth, '=');
  progressBar += std::string(remainingWidth, ' ');
  progressBar += "] " + std::to_string(floor(percent)) + "%";
  
  Rcpp::Rcout << progressBar;
}