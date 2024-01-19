#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "Circle.h"
#include "ClusterList.h"
#include "math.h"

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

// helper class for rcppConstructLineLinkDf
class LineLinkDataFrameFactory {
private:
    std::vector<ClusterList> clusterListVector;
    int numClusters;

    std::vector<std::unordered_map<std::string, int>> clusteredClonotypeIndex;
    int numClusteredClonotypes;

    std::vector<double> x1, x2, y1, y2;
    std::vector<int> cluster1, cluster2;

public:
    static Rcpp::DataFrame constructFrom(
        Rcpp::List clusterLists,
        Rcpp::List rawCloneSizes,
        Rcpp::List sharedClonotypeClusters,
        double extraSpacing
    ) {
        LineLinkDataFrameFactory llDfFactory (
            clusterLists, rawCloneSizes, sharedClonotypeClusters, extraSpacing
        );
        return llDfFactory.createOutputDataFrame();
    }

private:
    LineLinkDataFrameFactory(
        Rcpp::List clusterLists,
        Rcpp::List rawCloneSizes,
        Rcpp::List sharedClonotypeClusters,
        double extraSpacing
    ) {
        clusterListVector = createCppClusterListVector(clusterLists);
        clusteredClonotypeIndex = createClusteredClonotypeIndex(rawCloneSizes);

        std::vector<std::string> clonotypes (sharedClonotypeClusters.names());
        std::vector<std::vector<int>> clusterIndicies = Rcpp::as<std::vector<std::vector<int>>>(
            sharedClonotypeClusters
        );

        for (int i = 0; i < (int) clonotypes.size(); i++) {

            int currNumSharedClusters = clusterIndicies[i].size();
            std::vector<Circle> currCircles (currNumSharedClusters);
            std::vector<int> currOneIndexedClusterIndicies (currNumSharedClusters);

            for (int j = 0; j < currNumSharedClusters; j++) {
                ClusterList& currSharedCluster = clusterListVector[clusterIndicies[i][j]];
                currCircles[j] = currSharedCluster.getClonotypeCircle(clonotypes[i]);
                currOneIndexedClusterIndicies[j] = clusterIndicies[i][j];
            }

            addSharedCircleLinkInfo(currCircles, currOneIndexedClusterIndicies, extraSpacing);
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
            if (currentClonotypeTable.size() == 0) {continue;};

            addElementsToHashMap(
                outputIndex[i],
                Rcpp::as<std::vector<std::string>>(currentClonotypeTable.names()),
                Rcpp::as<std::vector<int>>(currentClonotypeTable)
            );
        }

        return outputIndex;
    }

    static void addElementsToHashMap(
        std::unordered_map<std::string, int>& hashMap,
        std::vector<std::string> keys,
        std::vector<int> values
    ) {
        for (int i = 0; i < (int) keys.size(); i++) {
            hashMap.emplace(keys[i], values[i]);
        }
    }

    // this is dependent on if the user wants to show every link
    void addSharedCircleLinkInfo(
        std::vector<Circle>& circles,
        std::vector<int>& currOneIndexedClusterIndicies,
        double extraSpacing
    ) {
        int numCircles = circles.size();
        for (int i = 0; i < numCircles - 1; i++) {
            for (int j = i + 1; j < numCircles; j++) {
                addTwoSharedCircleLinkInfo(circles[i], circles[j], extraSpacing);

                // push the indicies : not guaranteed atm to be in correct order
                cluster1.push_back(currOneIndexedClusterIndicies[i] + 1);
                cluster2.push_back(currOneIndexedClusterIndicies[j] + 1);
            }
        }
    }

    void addTwoSharedCircleLinkInfo(Circle& c1, Circle& c2, double extraSpacing) {

        TwoDRightPointingLine linkLine = TwoDRightPointingLineFactory::createCircleLinkingLineWithSpacing(
            c1, c2, extraSpacing
        );
        
        x1.push_back(linkLine.getLeftX());
        y1.push_back(linkLine.getLeftY());

        x2.push_back(linkLine.getRightX());
        y2.push_back(linkLine.getRightY());

    }

    Rcpp::DataFrame createOutputDataFrame() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("x1") = x1,
            Rcpp::Named("x2") = x2,
            Rcpp::Named("y1") = y1,
            Rcpp::Named("y2") = y2,
            Rcpp::Named("c1") = cluster1,
            Rcpp::Named("c2") = cluster2
        );
    }
};

// [[Rcpp::export]]
Rcpp::DataFrame rcppConstructLineLinkDf(
    Rcpp::List clusterLists,
    Rcpp::List rawCloneSizes,
    Rcpp::List sharedClonotypeClusters,
    double extraSpacing
) {
    return LineLinkDataFrameFactory::constructFrom(
        clusterLists, rawCloneSizes, sharedClonotypeClusters, extraSpacing
    );
}
