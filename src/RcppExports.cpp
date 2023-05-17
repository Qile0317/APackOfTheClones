// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// estimate_rad
double estimate_rad(std::vector<double> x_vals, std::vector<double> rad_vals, double centroid_x);
RcppExport SEXP _APackOfTheClones_estimate_rad(SEXP x_valsSEXP, SEXP rad_valsSEXP, SEXP centroid_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x_vals(x_valsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rad_vals(rad_valsSEXP);
    Rcpp::traits::input_parameter< double >::type centroid_x(centroid_xSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_rad(x_vals, rad_vals, centroid_x));
    return rcpp_result_gen;
END_RCPP
}
// get_transformed_clone_sizes
Rcpp::List get_transformed_clone_sizes(Rcpp::List sizelist, double clone_scale_factor, int num_clusters);
RcppExport SEXP _APackOfTheClones_get_transformed_clone_sizes(SEXP sizelistSEXP, SEXP clone_scale_factorSEXP, SEXP num_clustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type sizelist(sizelistSEXP);
    Rcpp::traits::input_parameter< double >::type clone_scale_factor(clone_scale_factorSEXP);
    Rcpp::traits::input_parameter< int >::type num_clusters(num_clustersSEXP);
    rcpp_result_gen = Rcpp::wrap(get_transformed_clone_sizes(sizelist, clone_scale_factor, num_clusters));
    return rcpp_result_gen;
END_RCPP
}
// get_average_vector
std::vector<double> get_average_vector(Rcpp::List vec_list);
RcppExport SEXP _APackOfTheClones_get_average_vector(SEXP vec_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type vec_list(vec_listSEXP);
    rcpp_result_gen = Rcpp::wrap(get_average_vector(vec_list));
    return rcpp_result_gen;
END_RCPP
}
// get_component_repulsion_vector
std::vector<double> get_component_repulsion_vector(Rcpp::List inp, int i, int j, double G);
RcppExport SEXP _APackOfTheClones_get_component_repulsion_vector(SEXP inpSEXP, SEXP iSEXP, SEXP jSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type inp(inpSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(get_component_repulsion_vector(inp, i, j, G));
    return rcpp_result_gen;
END_RCPP
}
// do_cluster_intersect
bool do_cluster_intersect(std::vector<double> Cn_centroid, double Cn_clRad, std::vector<double> Cm_centroid, double Cm_clRad, double thr);
RcppExport SEXP _APackOfTheClones_do_cluster_intersect(SEXP Cn_centroidSEXP, SEXP Cn_clRadSEXP, SEXP Cm_centroidSEXP, SEXP Cm_clRadSEXP, SEXP thrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type Cn_centroid(Cn_centroidSEXP);
    Rcpp::traits::input_parameter< double >::type Cn_clRad(Cn_clRadSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Cm_centroid(Cm_centroidSEXP);
    Rcpp::traits::input_parameter< double >::type Cm_clRad(Cm_clRadSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    rcpp_result_gen = Rcpp::wrap(do_cluster_intersect(Cn_centroid, Cn_clRad, Cm_centroid, Cm_clRad, thr));
    return rcpp_result_gen;
END_RCPP
}
// calculate_transformation_vectors
Rcpp::List calculate_transformation_vectors(Rcpp::List transformation_vectors, Rcpp::List overall_repulsion_vec, int num_clusters);
RcppExport SEXP _APackOfTheClones_calculate_transformation_vectors(SEXP transformation_vectorsSEXP, SEXP overall_repulsion_vecSEXP, SEXP num_clustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type transformation_vectors(transformation_vectorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type overall_repulsion_vec(overall_repulsion_vecSEXP);
    Rcpp::traits::input_parameter< int >::type num_clusters(num_clustersSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_transformation_vectors(transformation_vectors, overall_repulsion_vec, num_clusters));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_APackOfTheClones_estimate_rad", (DL_FUNC) &_APackOfTheClones_estimate_rad, 3},
    {"_APackOfTheClones_get_transformed_clone_sizes", (DL_FUNC) &_APackOfTheClones_get_transformed_clone_sizes, 3},
    {"_APackOfTheClones_get_average_vector", (DL_FUNC) &_APackOfTheClones_get_average_vector, 1},
    {"_APackOfTheClones_get_component_repulsion_vector", (DL_FUNC) &_APackOfTheClones_get_component_repulsion_vector, 4},
    {"_APackOfTheClones_do_cluster_intersect", (DL_FUNC) &_APackOfTheClones_do_cluster_intersect, 5},
    {"_APackOfTheClones_calculate_transformation_vectors", (DL_FUNC) &_APackOfTheClones_calculate_transformation_vectors, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_APackOfTheClones(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}