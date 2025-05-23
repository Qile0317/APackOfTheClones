#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "LineLinkDataFrameFactory.h"

// TODO make class for shared clones in C++

// helper for the remove_unique_clones function
// returns a length 2 list with the first element being a vector of clonotypes
// and the second element being a list of clusters
// [[Rcpp::export]]
Rcpp::List rcppRemoveUniqueClonesHelper(
    std::vector<std::string> clonotypes, std::vector<std::vector<int>> clusters
) {
    std::vector<std::string> filteredClonotypes;
    std::vector<std::vector<int>> filteredClusters;

    for (int i = 0; i < (int) clusters.size(); i++) {
        if ((int) clusters[i].size() <= 1) {
            continue;
        }
        filteredClonotypes.push_back(clonotypes[i]);
        filteredClusters.push_back(clusters[i]);
    }

    return Rcpp::List::create(filteredClonotypes, filteredClusters);
}

// [[Rcpp::export]]
Rcpp::DataFrame rcppConstructLineLinkDf(
    Rcpp::List clusterLists,
    Rcpp::List rawCloneSizes,
    Rcpp::List sharedClonotypeClusters,
    int oneIndexedSourceClusterIndex,
    double extraSpacing,
    bool showAllLinks
) {
    return LineLinkDataFrameFactory::constructFrom(
        clusterLists,
        rawCloneSizes,
        sharedClonotypeClusters,
        oneIndexedSourceClusterIndex,
        extraSpacing,
        showAllLinks
    );
}
