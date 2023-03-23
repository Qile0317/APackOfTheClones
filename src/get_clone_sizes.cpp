#include <Rcpp.h>
#include <RcppCommon.h>
#include <unordered_map>
#include <vector>
#include <string>

using namespace Rcpp;

typedef std::unordered_map<String, double> cluster_clonotypes_hash; // double has about 15-16 digits of precision so this should be fine
typedef std::vector<cluster_clonotypes_hash> cluster_clonotypes_hash_vec;
typedef std::vector<NumericVector> result_vector_list;

// convinience function to initialize a vector of empty hashmaps
cluster_clonotypes_hash_vec init_vec(int& num_clusters) {
  cluster_clonotypes_hash_vec output_vec;
  for (int i = 0; i < num_clusters; i++) {
    cluster_clonotypes_hash empty_map = {};
    output_vec.push_back(empty_map);
  }
  return output_vec;
}

// convinience function to initalize the result list.
// Apparently a C++ vector of R numeric vectors returns a list of numeric vectors in R
result_vector_list init_result_vec(int& num_clusters) {
  result_vector_list output_list;
  for (int i = 0; i < num_clusters; i++) {
    NumericVector empty_vec = NumericVector::create();
    output_list.push_back(empty_vec);
  }
  return output_list;
}

// convert vector of hashmaps into the R output list of numeric vectors
result_vector_list convert_to_list_output(cluster_clonotypes_hash_vec vector_of_hashmap_of_clonotype_count_per_cluster, int& num_clusters) {
  result_vector_list list_output = init_result_vec(num_clusters);
  for (int curr_cluster = 0; curr_cluster < num_clusters; curr_cluster++) {
    for (auto& key_val_pair: vector_of_hashmap_of_clonotype_count_per_cluster[curr_cluster]) {
      list_output[curr_cluster].push_back(key_val_pair.second);
    }
  }
  return list_output;
}

//' @export
// [[Rcpp::export]]
std::vector<NumericVector> get_clone_sizes_Cpp(StringVector barcodes, NumericVector clusters, StringVector clonotype_ids, int num_clusters, double scale_factor) {

  result_vector_list result_list_of_vector = init_result_vec(num_clusters);
  cluster_clonotypes_hash_vec vector_of_hashmap_of_clonotype_count_per_cluster = init_vec(num_clusters);
  int curr_cluster;

  // do the clonotype/cluster counting into a vector of hashmaps
  for (int i = 0; i < barcodes.length(); i++) {
    if (!StringVector::is_na(barcodes[i])) {
      String curr_clonotype = clonotype_ids[i];
      curr_cluster = clusters[i] - 1;

      if (vector_of_hashmap_of_clonotype_count_per_cluster[curr_cluster].count(curr_clonotype)) {
        vector_of_hashmap_of_clonotype_count_per_cluster[curr_cluster][curr_clonotype] += scale_factor;
      } else {
        vector_of_hashmap_of_clonotype_count_per_cluster[curr_cluster][curr_clonotype] = scale_factor;
      }
    }
  }
  return convert_to_list_output(vector_of_hashmap_of_clonotype_count_per_cluster, num_clusters);
}
