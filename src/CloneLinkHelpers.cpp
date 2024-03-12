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

// function to filter shared clones. outputs list(filtered_names, filtered_clusters)
// [[Rcpp::export]]
Rcpp::List rcppFilterSharedClonesByClusterHelper(
    std::vector<std::vector<int>> sharedClusters, // one-indexed!
    std::vector<std::string> clonotypes,
    std::vector<bool> includeCluster
) {
    std::vector<std::vector<int>> filteredSharedClusters;
    std::vector<std::string> filteredclonotypes;

    for (int i = 0; i < (int) sharedClusters.size(); i++) {
        for (int clusterIndex : sharedClusters[i]) {
            if (!includeCluster[clusterIndex - 1]) {continue;}
            filteredSharedClusters.push_back(sharedClusters[i]);
            filteredclonotypes.push_back(clonotypes[i]);
            break;
        }
    }

    return Rcpp::List::create(filteredclonotypes, filteredSharedClusters);
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
        Rcpp::List rawCloneSizes, // list of table objects, some may be empty
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
            if (currentClonotypeTable.size() == 0) {continue;}

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
        for (int i = 0; i < ((int) circles.size()) - 1; i++) {
            for (int j = i + 1; j < (int) circles.size(); j++) {
                
                // TODO for filtered ver, only add info if one of circles is from an origin cluster.

                TwoDLine linkLine = TwoDLineFactory::createCircleLinkingLineWithSpacing(
                    circles[i], circles[j], extraSpacing
                );

                x1.push_back(linkLine.getLeftX());
                y1.push_back(linkLine.getLeftY());

                x2.push_back(linkLine.getRightX());
                y2.push_back(linkLine.getRightY());

                int leftCircleIndex = linkLine.matchLeftCircleIndex(circles, i, j);
                int rightCircleIndex = i + j - leftCircleIndex;

                cluster1.push_back(currOneIndexedClusterIndicies[leftCircleIndex]);
                cluster2.push_back(currOneIndexedClusterIndicies[rightCircleIndex]);

            }
        }
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
