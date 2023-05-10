#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

inline double sqr(double n) {
  return n*n;
}

inline std::vector<double> neg_dir_vec(List c1, List c2) {
  std::vector<double> c1_centroid = c1[3], c2_centroid = c2[3];
  return {c1_centroid[0] - c2_centroid[0], c1_centroid[1] - c2_centroid[1]};
}

inline std::vector<double> polar_form(std::vector<double> v) {
  return {sqrt(sqr(v[0]) + sqr(v[1])), atan2(v[1], v[0])};
}

// [[Rcpp::export]]
std::vector<double> pdV(List c1, List c2) {
  return polar_form(neg_dir_vec(c1, c2));
}

// [[Rcpp::export]]
std::vector<double> comV(std::vector<double> polar_vec) {
  return {polar_vec[0] * cos(polar_vec[1]), polar_vec[0] * sin(polar_vec[1])};
}

// [[Rcpp::export]]
std::vector<double> get_average_vector(List vec_list) {
  std::vector<double> sum_vector = {0, 0};
  int n = vec_list.size();
  for (int i = 0; i < n; i++) {
    std::vector<double> curr_vec = vec_list[i];
    sum_vector[0] += curr_vec[0];
    sum_vector[1] += curr_vec[1];
  }
  sum_vector[0] /= n;
  sum_vector[1] /= n;
  return sum_vector;
}