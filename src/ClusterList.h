#pragma once

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "Circle.h"
#include "math.h"

class ClusterList {
protected: // variables

    std::unordered_map<std::string, int> clonotypeIndex;
    std::vector<Circle> circles;
    double centroidX;
    double centroidY;
    double clRad;
    int numClones;
    bool isEmpty;

public: // constructor

    ClusterList(const Rcpp::List rClusterList) {
       
        if (rClusterList.size() == 0) {
            isEmpty = true;
            return;
        }

        isEmpty = false;
        clRad = rClusterList["clRad"];
        centroidX = getRCentroidX(rClusterList);
        centroidY = getRCentroidY(rClusterList);
        clonotypeIndex = rCharactersToHashMap(rClusterList["clonotype"]);

        // get x,y,r and make into circles
        Rcpp::NumericVector x = rClusterList["x"], y = rClusterList["y"], r = rClusterList["rad"];
        numClones = x.size();

        for (int i = 0; i < numClones; i++) {
            circles.push_back(Circle(x[i], y[i], r[i]));
        }

    }

private:

    double getRCentroidX(Rcpp::List rClusterList) {
        Rcpp::NumericVector centroidVector = rClusterList["centroid"];
        return centroidVector[0];
    }

    double getRCentroidY(Rcpp::List rClusterList) {
        Rcpp::NumericVector centroidVector = rClusterList["centroid"];
        return centroidVector[1];
    }

    std::unordered_map<std::string, int> rCharactersToHashMap(Rcpp::CharacterVector v) {
        std::unordered_map<std::string, int> outputHashMap (v.size());
        for (int i = 0; i < (int) v.size(); i++) {
            outputHashMap.emplace(v[i], i);
        }
        return outputHashMap;
    }

public: // static methods

    static double toRadDecrease(double cloneScale, double radScale) {
        return cloneScale * (1 - radScale);
    }

public: // general methods

    bool isEmptyClusterList() {
        return isEmpty;
    }

    Circle getClonotypeCircle(const std::string clonotype) {
        return circles[clonotypeIndex[clonotype]];
    }

    ClusterList& rescaleClones(double newCloneScale, double prevCloneScale, double prevRadScale) {

        if (isEmpty) return *this;

        double scaleFactor = newCloneScale / prevCloneScale;
        double prevRadDecrease = toRadDecrease(prevCloneScale, prevRadScale);
        double newRadDecrease = toRadDecrease(newCloneScale, prevRadScale);

        for (int i = 0; i < (int) circles.size(); i++) {

            TwoDVector newOrigin = TwoDVector::createCircleOrigin(circles[i])
                .decreaseXYComponent(centroidX, centroidY)
                .scaleMagnitude(scaleFactor)
                .increaseXYComponent(centroidX, centroidY);

            circles[i]
                .setX(newOrigin.getX())
                .setY(newOrigin.getY())
                .increaseRad(prevRadDecrease)
                .scaleRad(scaleFactor)
                .decreaseRad(newRadDecrease);

        }

        return recalculateClusterRadius(newRadDecrease);
    }

private:

    ClusterList& recalculateClusterRadius(double radDecrease) {

        if (circles.size() == 1) {
            clRad = circles[0].rad() + radDecrease;
            return *this;
        }

        if (circles.size() == 2) {
            clRad = (0.5 * (circles[0].rad() + circles[1].rad()));
            return *this;
        }

        // imperfect cluster radius approximator - same as CirclePacker
        int xmaxIndex = 0;
        for (int i = 1; i < (int) circles.size(); i++) {
            if (circles[i].x() > circles[xmaxIndex].x()) {
                xmaxIndex = i;
            }
        }
        clRad = circles[xmaxIndex].x() + circles[xmaxIndex].rad() - centroidX;
        return *this;
    }

public: // getters

    Rcpp::List getRClusterList() {

        if (isEmpty) {
            return Rcpp::List::create();
        }

        return Rcpp::List::create(
            Rcpp::Named("x") = getXVector(),
            Rcpp::Named("y") = getYVector(),
            Rcpp::Named("rad") = getRadiiVector(),
            Rcpp::Named("centroid") = getRCentroidVector(),
            Rcpp::Named("clRad") = getApproximatedClusterRadius(),
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
        return Rcpp::NumericVector::create(centroidX, centroidY);
    }

    double getCentroidX() {
        return centroidX;
    }

    double getCentroidY() {
        return centroidY;
    }

    double getApproximatedClusterRadius() {
        return clRad;
    }

    std::vector<std::string> getClonotypeVector() {
        std::vector<std::string> clonotypes(numClones);
        for (auto keyValuePair : clonotypeIndex) {
            clonotypes[keyValuePair.second] = keyValuePair.first;
        }
        return clonotypes;
    }

};
