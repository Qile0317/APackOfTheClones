using Rcpp::_; // for overlap_check

#define all(x) (x).begin(), (x).end()
#define Llui long long unsigned int
#define Rvec Rcpp::Nullable<Rcpp::NumericVector>
#define _progbar(b,c) if(verbose) {progress_bar(b,c);}

inline double sqr(double n) {
    return n*n;
}

bool approx_equal(double a, double b, double epsilon = 5e-5) {
    return std::abs(a - b) <= epsilon;
}

// function to compare two numericvectors of the same length
bool elements_are_equal(
    Rcpp::NumericVector vec1, Rcpp::NumericVector vec2, double epsilon = 5e-5
) {
    for (int i = 0; i < vec1.size(); i++) {
        if (!approx_equal(vec1[i], vec2[i], epsilon)) {
            return false;
        }
    }
    return true;
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

// simple alias for an user interupt checker
inline void periodic_interupt_check(int num, int factor) {
    if ((num % factor) == 0) {
      Rcpp::checkUserInterrupt();
    }
}