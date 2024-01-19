#include <Rcpp.h>
#include <vector>

// Function to subset to only edge circles in C++
// [[Rcpp::export]]
std::vector<int> rcppGetEdgeCircleIndicies(Rcpp::DataFrame apotcPlotDataFrame) {

    // target indicies
    int xmin = 0, xmax = 0, ymin = 0, ymax = 0;

    // get relevant dataframe column references
    Rcpp::NumericVector x = apotcPlotDataFrame["x"];
    Rcpp::NumericVector y = apotcPlotDataFrame["y"];
    Rcpp::NumericVector r = apotcPlotDataFrame["r"];

    for (int i = 0; i < (int) x.size(); i++) {
        double xi = x[i], yi = y[i], ri = r[i];

        if (xi - ri < x[xmin] - r[xmin]) {
            xmin = i;
        }
        if (xi + ri > x[xmax] + r[xmax]) {
            xmax = i;
        }
        if (yi - ri < y[ymin] - r[ymin]) {
            ymin = i;
        }
        if (yi + ri > y[ymax] + r[ymax]) {
            ymax = i;
        }
    }

    return {xmin + 1, xmax + 1, ymin + 1, ymax + 1};
}
