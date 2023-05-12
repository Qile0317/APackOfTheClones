#include <Rcpp.h>
#include <vector>
#include <cmath>

using Rvec = Rcpp::Nullable<Rcpp::NumericVector>;

inline double sqr(double n) {
  return n*n;
}

std::vector<double> neg_dir_vec(Rcpp::List c1, Rcpp::List c2) {
  std::vector<double> c1_centroid = c1[3], c2_centroid = c2[3];
  return {c1_centroid[0] - c2_centroid[0], c1_centroid[1] - c2_centroid[1]};
}

std::vector<double> polar_form(std::vector<double> v) {
  return {sqrt(sqr(v[0]) + sqr(v[1])), atan2(v[1], v[0])};
}

std::vector<double> neg_polar_dist_vec(Rcpp::List c1, Rcpp::List c2) {
  return polar_form(neg_dir_vec(c1, c2));
}

std::vector<double> component_form(std::vector<double> polar_vec) {
  return {polar_vec[0] * cos(polar_vec[1]), polar_vec[0] * sin(polar_vec[1])};
}

// calculate avg vector of an R Rcpp::List of 2D vectors, skipping vectors of 0, 0
// [[Rcpp::export]]
std::vector<double> get_average_vector(Rcpp::List vec_list) {
  std::vector<double> sum_vector = {0, 0};
  std::vector<double> blank = {0, 0};
  int num_non_zero_vectors = 0;
  for (int i = 0; i < vec_list.size(); i++) {
    std::vector<double> curr_vec = vec_list[i];
    if (curr_vec != blank) {
      num_non_zero_vectors++;
      sum_vector[0] += curr_vec[0];
      sum_vector[1] += curr_vec[1];
    }
  }
  if (num_non_zero_vectors) {
    sum_vector[0] /= num_non_zero_vectors;
    sum_vector[1] /= num_non_zero_vectors;
  }
  return sum_vector;
}

// compute component form of repulsion vector between two clusters in the
// clusterRcpp::List `inp`, assuming do_proceed(inp,i,j,thr) == TRUE
// uses a modified version of coulombs law except its addition in the numerator
// to compute the magnitude of the repulsion vector / 2 with the same direction
// [[Rcpp::export]]
std::vector<double> get_component_repulsion_vector(
  Rcpp::List inp, int i, int j, double G
) {
  Rcpp::List c1 = inp[i-1], c2 = inp[j-1];
  double c1_rad = c1[4], c2_rad = c2[4];
  std::vector<double> neg_polar_vec = neg_polar_dist_vec(c1, c2);
  std::vector<double> polar_repulsion_vec = {0, neg_polar_vec[1]};
  polar_repulsion_vec[0] = 0.5 * G * (c1_rad + c2_rad) / sqr(neg_polar_vec[0]);
  return component_form(polar_repulsion_vec);
} 

// TODO: initialize_direction_vectors, initialize_list_of_transformation_vectors

// Check if 2 cluster lists overlap, with a threshold, give their centroids and
// radii. thr should be the amount of acceptable overlap
// [[Rcpp::export]]
bool do_cluster_intersect(
  std::vector<double> Cn_centroid, double Cn_clRad,
  std::vector<double> Cm_centroid, double Cm_clRad,
  double thr
) {
  double x_dif = sqr(Cn_centroid[0] - Cm_centroid[0]);
  double y_dif = sqr(Cn_centroid[1] - Cm_centroid[1]);
  return (sqrt(x_dif + y_dif) + thr) < (Cn_clRad + Cm_clRad);
}

bool do_proceed(Rcpp::List inp, int i, int j, double thr) {
  if (i != j) {
    Rcpp::List Cn = inp[i], Cm = inp[j];
    if (Cn.size()) {
      if (Cm.size()) {
        std::vector<double> Cn_centroid = Cn[3], Cm_centroid = Cm[3];
        double Cn_rad = Cn[4], Cm_rad = Cm[4];
        return do_cluster_intersect(
          Cn_centroid, Cn_rad, Cm_centroid, Cm_rad, thr
        );
      }
    }
  }
  return false;
}

/* O(N^2) operation to calculate all repulsion vectors for each cluster in list

Rcpp::List calculate_repulsion_vectors(
  Rcpp::List overall_repulsion_vec,
  Rcpp::List inp,
  int num_clusters,
  double G = 1,
  double thr = 0
) {
  for (int i = 0; i < num_clusters; i++) {
    
    for (int j = 0; j < num_clusters; j++) {
      
      if (do_proceed(inp, i, j, thr)) {
        
      }
    }
  }
}
*/

// function to just get the average vectors from a list of list of repulsion
// vectors within 1 iteration
// [[Rcpp::export]]
Rcpp::List calculate_transformation_vectors(
    Rcpp::List transformation_vectors,
    Rcpp::List overall_repulsion_vec,
    int num_clusters
) {
  bool contains_nonzero_vector = false;
  std::vector<double> blank = {0, 0};
  for (int i = 0; i < num_clusters; i++) {
    Rcpp::List curr_rep_vecs = overall_repulsion_vec[i];
    transformation_vectors[i] = get_average_vector(curr_rep_vecs);
    std::vector<double> tv = transformation_vectors[i];
    if (tv != blank) {
      contains_nonzero_vector = true;
    }
  }
  if (contains_nonzero_vector) {
    return transformation_vectors;
  }
  return Rcpp::List::create();
}