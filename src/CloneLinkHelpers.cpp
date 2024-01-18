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

// [[Rcpp::export]]
Rcpp::DataFrame rcppConstructLineLinkDf(
    Rcpp::List clusterLists,
    Rcpp::List rawCloneSizes,
    std::vector<std::vector<int>> sharedClones
) {
    return LineLinkDataFrameFactory::constructFrom(
        clusterLists, rawCloneSizes, sharedClones
    )
}

class LineLinkDataFrameFactory {
private:
    std::vector<ClusterList> clusterListVector;
    int numClusters;

    std::vector<&std::unordered_map<std::string, int>> clusteredClonotypeIndex;
    int numClusteredClonotypes;

    std::vector<double> x1, x2, y1, y2, r1, r2;
    int numColumns;

public:
    static Rcpp::DataFrame constructFrom(
        Rcpp::List clusterLists,
        Rcpp::List rawCloneSizes,
        std::vector<std::vector<int>> sharedClones
    ) {
        LineLinkDataFrameFactory factory ();
        return factory.createOutputDataFrame();
    }

private:
    LineLinkDataFrameFactory(
        Rcpp::List clusterLists,
        Rcpp::List rawCloneSizes,
        Rcpp::List sharedClonotypeClusters,
    ) {
        clusterListVector = createCppclusterListVector(clusterLists);
        clusteredClonotypeIndex = writeClusteredClonotypeIndex(rawCloneSizes);

        std::vector<std::string> clonotypes (sharedClonotypeClusters.names());
        std::vector<std::vector<int>> clusterIndicies = Rcpp::as<std::vector<std::vector<int>>>(
            sharedClonotypeClusters
        );

        for (int i = 0; i < (int) clonotypes.size(); i++) {
            std::vector<Circle> currCircles (currNumSharedClusters);
            for (int j = 0; j < (int) clonotypes[i].size(); j++) {
                currCircles.push_back(
                    clusterListVector[clusterIndicies[j]].getClonotypeCircle(clonotypes[i])
                )
            }
            addSharedCircleLinkInfo(currCircles);
        }
    }

    std::vector<ClusterList> createCppClusterListVector(Rcpp::List clusterLists) {
        numClusters = clusterLists.size(); // assume non-zero
        std::vector<ClusterList> outputClusterListVector (numClusters);

        for (int i = 0; i < numClusters; i++) {
            outputClusterListVector[i] = ClusterList(clusterLists[i]);
        }

        outputClusterListVector.resize(numClusters);
        return outputClusterListVector;
    }

    std::vector<std::unordered_map<std::string, int>> createClusteredClonotypeIndex(
        Rcpp::List rawCloneSizes
    ) {

        numClusters = rawCloneSizes.size();
        std::vector<std::unordered_map<std::string, int>> outputIndex (
            numClusters, std::unordered_map<std::string, int>() // uncertain if this constructor works
        );

        for (int i = 0; i < numClusters; i++) {

            Rcpp::NumericVector currentClonotypeTable = rawCloneSizes[i];

            if (currentClonotypeTable.size() == 0) {
                continue;
            };

            addElementsToHashMap(
                outputIndex[i],
                Rcpp::as<std::vector<std::string>>(v.names()),
                Rcpp::as<std::vector<int>>(v)
            )
        }

        return outputIndex;
    }

    static void addElementsToHashMap(
        std::unordered_map<std::string, int>& hashMap,
        std::vector<std::string> keys,
        std::vector<std::int> values
    ) {
        for (int i = 0; i < (int) keys.size(); i++) {
            hashMap.emplace(keys[i], values[i]);
        }
    }

    // this is dependent on if the user wants to show every link
    void addSharedCircleLinkInfo(std::vector<Circle> circles) {
        for (int i = 0; i < (int) circles.size() - 1; i++) {
            for (int j = i + 1; j < (int) circles.size(); j++) {
                addTwoSharedCircleLinkInfo(circles[i], circles[j]);
            }
        }
    }

    void addTwoSharedCircleLinkInfo(Circle& c1, Circle& c2) {
        x1.push_back(c1.x());
        y1.push_back(c1.y());
        r1.push_back(c1.rad());

        x2.push_back(c2.x());
        y2.push_back(c2.y());
        r2.push_back(c2.rad());
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
