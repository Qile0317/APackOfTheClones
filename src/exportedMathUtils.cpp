#include <Rcpp.h>
#include <vector>
#include <string>
#include "math.h"

// assumes x > 1
// [[Rcpp::export]]
std::vector<std::vector<int>> rcppGetUniquePairsUpTo(int x, bool oneIndexed) {

    std::vector<std::vector<int>> uniquePairList (x * (x - 1) / 2);

    int index = 0;
    for (int i = oneIndexed; i < x - 1 + oneIndexed; i++)  {
        for (int j = i + 1; j < x + oneIndexed; j++) {
            uniquePairList[index] = {i, j};
            index++;
        }
    }

    return uniquePairList;
}
