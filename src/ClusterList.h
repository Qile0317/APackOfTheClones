#pragma once

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "Circle.h"

class ClusterList {
protected: // variables

    std::vector<Circle> circles;
    std::pair<double, double> centroid;
    double clRad;
    int numClones;

    std::unordered_map<std::string, int> clonotypeIndex;
    bool isEmpty;

public: // constructors

    ClusterList() {
        isEmpty = true;
    }

    ClusterList(const Rcpp::List rClusterList) {
        
        if (rClusterList.size() == 0) {
            isEmpty = true;
            return;
        }

        isEmpty = false;
        clRad = rClusterList["clRad"];
        centroid = getRCentroid(rClusterList);

        // get x,y,r and make into circles
        Rcpp::NumericVector x = rClusterList["x"], y = rClusterList["y"], r = rClusterList["rad"];
        numClones = x.size();

        for (int i = 0; i < numClones; i++) {
            circles.push_back(Circle(x[i], y[i], r[i]));
        }

        clonotypeIndex = std::unordered_map<std::string, int>();
        std::vector<std::string> clonotypes = rClusterList["clonotype"];
        for (int i = 0; i < numClones; i++) {
            clonotypeIndex.emplace(clonotypes[i], i);
        }
    }

private:

    std::pair<double, double> getRCentroid(Rcpp::List rClusterList) {
        Rcpp::NumericVector centroidVector = rClusterList["centroid"];
        return std::make_pair(centroidVector[0], centroidVector[1]);
    }

public:

    bool isEmptyClusterList() {
        return isEmpty;
    }

    Circle getClonotypeCircle(const std::string clonotype) {
        return circles[clonotypeIndex[clonotype]];
    }

    // getters

    Rcpp::List getRClusterList() {
        if (isEmptyClusterList()) {
            return Rcpp::List::create();
        }

        return Rcpp::List::create(
            Rcpp::Named("x") = getXVector(),
            Rcpp::Named("y") = getYVector(),
            Rcpp::Named("rad") = getRadiiVector(),
            Rcpp::Named("clRad") = getApproximatedClusterRadius(),
            Rcpp::Named("centroid") = getRCentroidVector(),
            Rcpp::Named("clonotype") = getClonotypeVector()
        );
    }

    std::vector<double> getXVector() {
        std::vector<double> x(numClones);
        for (int i = 0; i < numClones; i++) {
            x[i] = circles[i].x();
        }
        return x;
    }

    std::vector<double> getYVector() {
        std::vector<double> y(numClones);
        for (int i = 0; i < numClones; i++) {
            y[i] = circles[i].y();
        }
        return y;
    }

    std::vector<double> getRadiiVector() {
        std::vector<double> r(numClones);
        for (int i = 0; i < numClones; i++) {
            r[i] = circles[i].rad();
        }
        return r;
    }

    Rcpp::NumericVector getRCentroidVector() {
        return Rcpp::NumericVector::create(centroid.first, centroid.second);
    }

    double getApproximatedClusterRadius() {
        return clRad;
    }

    std::vector<std::string> getClonotypeVector() {
        std::vector<std::string> clonotypes(numClones);
        for (auto& [clonotype, index]: clonotypeIndex) {
            clonotypes[index] = clonotype;
        }
        return clonotypes;
    }
};
