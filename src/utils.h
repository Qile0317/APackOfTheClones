#include <Rcpp.h>
#include <string>

namespace utils {
    // functions for testing
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

    // [[Rcpp::export]]
    bool has_repeats(const Rcpp::NumericVector vec1, const Rcpp::NumericVector vec2) {
        std::unordered_set<double> s;

        for (int i = 0; i < vec1.size(); i++) {
            if (s.count(vec1[i]) > 0) {
                return true;
            }
            s.insert(vec1[i]);
        }

        for (int i = 0; i < vec2.size(); i++) {
            if (s.count(vec2[i]) > 0) {
                return true;
            }
            s.insert(vec2[i]);
        }
        return false;
    }

    // check if two string vectors has common elements or if the second vector has
    // repeated elements
    // [[Rcpp::export]]
    bool has_common_strs(
        const std::vector<std::string>& vec1, const std::vector<std::string>& vec2
    ) {
        // Create a hash set to store strings from vec1
        std::unordered_set<std::string> hashSet;
        for (const std::string& str : vec1) {
            hashSet.insert(str);
        }

        // Check if any string from vec2 is present in the hash set
        for (const std::string& str : vec2) {
            if (hashSet.count(str) > 0) {
                return true;
            }
            hashSet.insert(str);
        }
        return false;
    }
}
