#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

inline double sqr(double n) {
  return n*n;
}

std::vector<double> neg_dir_vec(List c1, List c2) {
  std::vector<double> c1_centroid = c1[3], c2_centroid = c2[3];
  return {c1_centroid[0] - c2_centroid[0], c1_centroid[1] - c2_centroid[1]};
}

std::vector<double> polar_form(std::vector<double> v) {
  return {sqrt(sqr(v[0]) + sqr(v[1])), atan2(v[1], v[0])};
}

std::vector<double> neg_polar_dist_vec(List c1, List c2) {
  return polar_form(neg_dir_vec(c1, c2));
}

std::vector<double> component_form(std::vector<double> polar_vec) {
  return {polar_vec[0] * cos(polar_vec[1]), polar_vec[0] * sin(polar_vec[1])};
}

// avg vector of a list of 2D vectors
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

// compute component form of repulsion vector between two clusters in the
// clusterlist `inp`, assuming do_proceed(inp,i,j,thr) == TRUE
// uses a modified version of coulombs law except its addition in the numerator
// to compute the magnitude of the repulsion vector / 2 with the same direction
// [[Rcpp::export]]
std::vector<double> get_component_repulsion_vector(
  List inp, int i, int j, double G
) {
  List c1 = inp[i-1], c2 = inp[j-1];
  double c1_rad = c1[4], c2_rad = c2[4];
  std::vector<double> neg_polar_vec = neg_polar_dist_vec(c1, c2);
  std::vector<double> polar_repulsion_vec = {0, neg_polar_vec[1]};
  polar_repulsion_vec[0] = 0.5*G*(c1_rad+c2_rad) / sqr(neg_polar_vec[0]);
  return component_form(polar_repulsion_vec);
}