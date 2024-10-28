#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>

// [[Rcpp::export]]
Rcpp::NumericMatrix rcppFilterReductionCoords(
    std::vector<std::string>& seuratBarcodes,
    Rcpp::NumericMatrix reductionCoords
) {
    int numBarcodes = (int) seuratBarcodes.size();
    Rcpp::NumericMatrix output (numBarcodes, 2);

    std::unordered_set<std::string> seuratBarcodeSet (
        seuratBarcodes.begin(), seuratBarcodes.end()
    );

    std::vector<std::string> reductionBarcodes = Rcpp::as<std::vector<std::string>>(
        Rcpp::rownames(reductionCoords)
    );
   
    int outputRowIndex = 0;
    for (int i = 0; i < (int) reductionBarcodes.size(); i++) {
        if (seuratBarcodeSet.find(reductionBarcodes[i]) == seuratBarcodeSet.end()) {
            continue;
        }
        output (outputRowIndex, 0) = reductionCoords (i, 0);
        output (outputRowIndex, 1) = reductionCoords (i, 1);
        outputRowIndex++;
    }

    return output;
}
