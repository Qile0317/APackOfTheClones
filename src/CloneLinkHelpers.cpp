#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "Circle.h"
#include "ClusterList.h"
#include "math.h"

// TODO make class for shared clones in C++

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

    std::vector<double> x1, x2, y1, y2;
    std::vector<int> cluster1, cluster2;

public:
    static Rcpp::DataFrame constructFrom(
        Rcpp::List clusterLists,
        Rcpp::List rawCloneSizes,
        Rcpp::List sharedClonotypeClusters,
        int oneIndexedSourceClusterIndex,
        double extraSpacing,
        bool showAllLinks
    ) {

        LineLinkDataFrameFactory llDfFactory (
            clusterLists,
            rawCloneSizes,
            sharedClonotypeClusters,
            oneIndexedSourceClusterIndex,
            extraSpacing,
            showAllLinks
        );

        return llDfFactory.createOutputDataFrame();
    }

private:
    // constructor that just creates the dataframe immediately
    LineLinkDataFrameFactory(
        Rcpp::List clusterLists, // list of potentially empty clusterlists
        Rcpp::List rawCloneSizes, // list of table objects, some may be empty
        Rcpp::List sharedClonotypeClusters, // output from getSharedClones
        int oneIndexedSourceClusterIndex,
        double extraSpacing,
        bool showAllLinks
    ) {
        clusterListVector = createCppClusterListVector(clusterLists);
        clusteredClonotypeIndex = createClusteredClonotypeIndex(rawCloneSizes);

        std::vector<std::string> clonotypes = sharedClonotypeClusters.names();
        std::vector<std::vector<int>> clusterIndices = getZeroIndexedClusterIndices(
            sharedClonotypeClusters
        );

        for (int i = 0; i < (int) clonotypes.size(); i++) {

            std::vector<Circle> currCircles;
            std::vector<int> currOneIndexedClusterIndices;

            for (int clusterIndex : clusterIndices[i]) {
                ClusterList& currSharedCluster = clusterListVector[clusterIndex];
                currCircles.push_back(currSharedCluster.getClonotypeCircle(clonotypes[i]));
                currOneIndexedClusterIndices.push_back(clusterIndex);
            }

            addSharedCircleLinkInfo(
                currCircles,
                currOneIndexedClusterIndices,
                extraSpacing,
                oneIndexedSourceClusterIndex,
                showAllLinks
            );
        }
    }

    std::vector<ClusterList> createCppClusterListVector(Rcpp::List clusterLists) {

        numClusters = clusterLists.size();
        std::vector<ClusterList> outputClusterListVector;
        outputClusterListVector.reserve(numClusters);

        for (int i = 0; i < numClusters; i++) {
            Rcpp::List currRClusterList = clusterLists[i];
            outputClusterListVector.push_back(ClusterList(currRClusterList));
        }
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

    std::vector<std::vector<int>> getZeroIndexedClusterIndices(
        Rcpp::List sharedClonotypeClusters
    ) {
        std::vector<std::vector<int>> output = Rcpp::as<std::vector<std::vector<int>>>(
            sharedClonotypeClusters
        );

        for (int i = 0; i < (int) output.size(); i++) {
            for (int j = 0; j < (int) output[i].size(); j++) {
                output[i][j] -= 1;
            }
        }

        return output;
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
        std::vector<int>& currOneIndexedClusterIndices,
        double extraSpacing,
        int oneIndexedSourceClusterIndex,
        bool showAllLinks // TODO: this should instead be a "strategy enum" that determines what i and j are used (thus how links are distributed)
    ) {

        auto processCirclePair = [&](int i, int j) {
            TwoDLine linkLine = TwoDLineFactory::createCircleLinkingLineWithSpacing(
                circles[i], circles[j], extraSpacing
            );

            int leftCircleIndex = linkLine.matchLeftCircleIndex(circles, i, j);
            int rightCircleIndex = linkLine.matchRightCircleIndex(circles, i, j);

            int leftCircleClusterIndex = currOneIndexedClusterIndices[leftCircleIndex];
            int rightCircleClusterIndex = currOneIndexedClusterIndices[rightCircleIndex];

            if (oneIndexedSourceClusterIndex != -1) {
                if (leftCircleClusterIndex + 1 != oneIndexedSourceClusterIndex
                    && rightCircleClusterIndex + 1 != oneIndexedSourceClusterIndex) {
                    return;
                }
            }

            cluster1.push_back(leftCircleClusterIndex);
            cluster2.push_back(rightCircleClusterIndex);
            x1.push_back(linkLine.getLeftX());
            y1.push_back(linkLine.getLeftY());
            x2.push_back(linkLine.getRightX());
            y2.push_back(linkLine.getRightY());
        };
        
        const int n = static_cast<int>(circles.size());
        for (int i = 0; i < n; i++) {
            for (int j = showAllLinks ? 0 : i + 1; j < n; j++) {
                if (i != j) {
                    processCirclePair(i, j);
                }
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
