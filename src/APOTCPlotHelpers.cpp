#include <Rcpp.h>
#include <vector>

// [[Rcpp::export]]
std::vector<int> rcppGetEdgeCircleindices(Rcpp::DataFrame apotcPlotDataFrame) {

    // target indices
    int xmin = 0, xmax = 0, ymin = 0, ymax = 0;

    // get relevant dataframe column references
    Rcpp::NumericVector x = apotcPlotDataFrame["x"];
    Rcpp::NumericVector y = apotcPlotDataFrame["y"];
    Rcpp::NumericVector r = apotcPlotDataFrame["r"];

    for (int i = 0; i < (int) x.size(); i++) {
        
        if (x[i] - r[i] < x[xmin] - r[xmin]) {
            xmin = i;
        }

        if (x[i] + r[i] > x[xmax] + r[xmax]) {
            xmax = i;
        }

        if (y[i] - r[i] < y[ymin] - r[ymin]) {
            ymin = i;
        }

        if (y[i] + r[i] > y[ymax] + r[ymax]) {
            ymax = i;
        }
    }

    return {xmin + 1, xmax + 1, ymin + 1, ymax + 1};
}
