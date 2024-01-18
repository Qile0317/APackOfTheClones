#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "Circle.h"

class ClusterList {
private:
    std::vector<Circle> circles;
    std::pair<double, double> centroid;
    double clRad;
    int numClones;

    // important
    std::unordered_map<std::string, int> clonotypeIndex;
    bool isEmpty;

public:
    ClusterList(const Rcpp::List rClusterList) {
        if (rClusterList.size() == 0) {
            isEmpty = true;
            return;
        }

        isEmpty = false;
        clRad = rClusterList["clRad"];
        centroid = getRCentroid(rClusterList);

        // get x,y,r and make into circles
        NumericVector x = rClusterList["x"], y = rClusterList["y"], r = rClusterList["rad"];
        numClones = x.size();

        for (int i = 0; i < numClones; i++) {
            circles.push_back(Circle(x[i], y[i], r[i]));
        }

        std::vector<std::string> clonotypes = rClusterList["clonotype"];
        for (int i = 0; i < numClones; i++) {
            clonotypeIndex[Rcpp::as<std::string>(clonotypes[i])] = i;
        }
    }

private:
    std::pair<double, double> getRCentroid(Rcpp::List rClusterList) {
        NumericVector centroidVector = rClusterList["centroid"];
        return std::make_pair(centroidVector[0], centroidVector[1]);
    }

public:

    bool isEmptyClusterList() {
        return isEmpty;
    }

    Circle getClonotypeCircle(std::string clonotype) {
        return getCircleAt([clonotypeIndex[clonotype]]);
    }

    Circle getCircleAt(int index) {
        return circles[index];
    }

    // getters
    double getEstimatedClusterRad() {
        return clRad;
    }

    std::pair<double, double> getCentroid() {
        return centroid;
    }

    int getNumClones() {
        return numClones;
    }
};
