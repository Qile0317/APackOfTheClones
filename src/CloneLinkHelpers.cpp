#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "Circle.h"
#include "ClusterList.h"

// helper for the remove_unique_clones function
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

    // ideally, return a named list with filteredClonotypes as names and filteredClusters as elements
    return Rcpp::List::create(filteredClonotypes, filteredClusters);
}

class LineLinkDataFrameFactory {
private:
    std::vector<&std::unordered_map<std::string, int>> clusteredClonotypeIndex;
    int numClusters;

    std::vector<double> x1, x2, y1, y2, r1, r2;
    int numColumns;

public:
    static Rcpp::DataFrame construct(
        Rcpp::List rawCloneSizes, // list of table objects
        std::vector<std::string> clonotypes,
        std::vector<std::vector<int>> sharedClones,
        std::string linkMode
    ) {
        LineLinkDataFrameFactory factory ();
        return factory.createOutputDataFrame();
    }

private:
    LineLinkDataFrameFactory(
        Rcpp::List clusterLists, // very messy, probably should have a clusterlist C++ wrapper
        Rcpp::List rawCloneSizes,
        Rcpp::List sharedClonotypeClusters,
    ) {
        writeClusteredClonotypeIndex(rawCloneSizes);

        std::vector<std::string> clonotypes (sharedClonotypeClusters.names());
        std::vector<std::vector<int>> clusterIndicies = Rcpp::as<std::vector<std::vector<int>>>(
            sharedClonotypeClusters
        );

        for (int i = 0; i < (int) clonotypes.size(); i++) {
            
        }
    }

    void writeClusteredClonotypeIndex(Rcpp::List rawCloneSizes) {

        numClusters = rawCloneSizes.size();

        for (int i = 0; i < numClusters; i++) {

            Rcpp::NumericVector currentClonotypeTable = rawCloneSizes[i];
            clusteredClonotypeIndex.push_back({});

            if (currentClonotypeTable.size() == 0) {
                continue;
            };

            addElementsToHashMap(
                clusteredClonotypeIndex[i],
                Rcpp::as<std::vector<std::string>>(v.names()),
                Rcpp::as<std::vector<int>>(v)
            )
            
        }
    }

    static void addElementsToHashMap(
        &std::unordered_map<std::string, int> hashMap,
        std::vector<std::string> keys,
        std::vector<std::int> values
    ) {
        for (int i = 0; i < (int) keys.size(); i++) {
            hashMap[keys[i]] = values[i];
        }
    }

    Rcpp::DataFrame createOutputDataFrame() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("x1") = x1,
            Rcpp::Named("x2") = x2,
            Rcpp::Named("y1") = y1,
            Rcpp::Named("y2") = y2,
            Rcpp::Named("r1") = r1,
            Rcpp::Named("r2") = r2
        )
    }
};

Rcpp::DataFrame rcppComputeLineLinkDf() {
    
}

// this works so tables can be caseted to numericvector :D
// [[Rcpp::export]]
void lol(Rcpp::NumericVector v) {
    std::vector<std::string> vn = v.names();
    for (int i = 0; i < v.size(); ++i)
        Rcpp::Rcout << vn[i] << " " << v[i] << "\n";
}