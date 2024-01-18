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
    int numClones

    // important
    std::unordered_map<std::string, int> clonotypeIndex;
    bool isEmpty;

public:
    ClusterList(Rcpp::List rClusterList) {
        if (rClusterList.size() == 0) {
            isEmpty = true;
            return;
        }

        isEmpty = false;
        clRad = rClusterList["clRad"];
        centroid = std::make_pair(rClusterList["centroid"][0], rClusterList["centroid"][1]);

        // unsure if compiler will optimize the accesses
        int numClones = rClusterList["x"].size();
        for (int i = 0; i < numClones; i++) {
            circles.push_back(
                Circle(rClusterList["x"][i], rClusterList["y"][i], rClusterList["rad"][i])
            );
        }

        std::vector<std::string> clonotypes = rClusterList["clonotype"];
        for (int i = 0; i < numClones; i++) {
            clonotypeIndex[Rcpp::as<std::string>(clonotypes[i])] = i;
        }
    }

    bool isEmpty() {
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
