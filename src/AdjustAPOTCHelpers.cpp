#include <Rcpp.h>
#include <vector>

#include "ClusterList.h"

// [[Rcpp::export]]
Rcpp::List rcppRescaleClones(
    Rcpp::List rClusterlist,
    double newCloneScale,
    double prevCloneScale,
    double prevRadScale
) {
    return ClusterList(rClusterlist)
        .rescaleClones(newCloneScale, prevCloneScale, prevRadScale)
        .getRClusterList();
}
