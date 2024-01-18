#pragma once

#include "Circle.h"

#include <cmath>

inline double sqr(double x) {
    return x * x;
}

class TwoDVector {
protected:
    double magnitude;
    double direction;

public:
    TwoDVector(double mag, double dir) : magnitude(mag), direction(dir) {}
    TwoDVector(TwoDVector& v) : magnitude(v.getMagnitude()), direction(v.getDirection()) {}

    // Overload +
    TwoDVector operator+(const TwoDVector& other) const {
        double dx = getX() + other.getX(), dy = getY() + other.getY();
        return TwoDVector(std::sqrt(sqr(dx) + sqr(dy)), std::atan2(dy, dx));
    }

    // Overload +=
    TwoDVector& operator+=(const TwoDVector& other) {
        *this = *this + other;
        return *this;
    }

    // Overload -
    TwoDVector operator-(const TwoDVector& other) const {
        double dx = getX() - other.getX(), dy = getY() - other.getY();
        return TwoDVector(std::sqrt(sqr(dx) + sqr(dy)), std::atan2(dy, dx));
    }

    double getMagnitude() const {
        return magnitude;
    }

    double getDirection() const {
        return direction;
    }

    double getX() const {
        return magnitude * std::cos(direction);
    }

    double getY() const {
        return magnitude * std::sin(direction);
    }

    TwoDVector unitVector() const {
        return TwoDVector(1, direction);
    }

    void setMagnitude(double val) {
        magnitude = val;
    }

    void decreaseMagnitude(double val) {
        magnitude -= val;
    }
};

class TwoDVectorFactory {
public:
    static TwoDVector createFromXY(double x, double y) {
        double mag = std::sqrt(sqr(x) + sqr(y));
        double dir = std::atan2(y, x);
        return TwoDVector(mag, dir);
    }

    static TwoDVector createCircleOrigin(Circle& circle) {
        return createFromXY(circle.x(), circle.y());
    }
};

class TwoDRightPointingLine {
protected:
    TwoDVector originCoordinate;
    TwoDVector lineVector; // points right of origin!

public:
    TwoDRightPointingLine(TwoDVector o, TwoDVector l) : originCoordinate(o), lineVector(l) {}

    double getLeftX() {
        return getLeftPoint().getX();
    }

    double getLeftY() {
        return getLeftPoint().getY();
    }

    TwoDVector getLeftPoint() {
        return originCoordinate;
    }

    double getRightX() {
        return getRightPoint().getX();
    }

    double getRightY() {
        return getRightPoint().getY();
    }

    TwoDVector getRightPoint() {
        return originCoordinate + lineVector;
    }

    // assumes their origins correspond to line origin and endpoints!
    void shortenByTwoCirclesRadiiAndExtraSpacing(
        Circle& c1, Circle& c2, double extraSpacing
    ) {
        Circle& leftCircle = c1, rightCircle = c2;
        if (c2.isLeftOf(c1)) {
            std::swap(leftCircle, rightCircle);
        }

        shortenByTwoCirclesRadii(
            leftCircle.rad() + extraSpacing, rightCircle.rad() + extraSpacing
        );
    }

    void shortenByTwoCirclesRadii(double leftRadius, double rightRadius) {
        shortenLeftByCircleRadius(leftRadius);
        shortenRightByCircleRadius(rightRadius);
    }

    void shortenLeftByCircleRadius(double val) {
        decreaseMagnitude(val);

        TwoDVector temp (lineVector);
        temp.setMagnitude(val);
        translate(temp);
    }

    void shortenRightByCircleRadius(double val) {
        decreaseMagnitude(val);
    }

    void decreaseMagnitude(double val) {
        lineVector.decreaseMagnitude(val);
    }

    void translate(TwoDVector v) {
        originCoordinate += v;
    }
};

class TwoDRightPointingLineFactory {
public:
    static TwoDRightPointingLine createCircleLinkingLineWithSpacing(
        Circle& c1, Circle& c2, double extraSpacing
    ) {
        Circle& leftCircle = c1, rightCircle = c2;
        if (c2.isLeftOf(c1)) {
            std::swap(leftCircle, rightCircle);
        }

        TwoDVector leftOrigin = TwoDVectorFactory::createCircleOrigin(leftCircle);

        TwoDVector rightOrigin = TwoDVectorFactory::createCircleOrigin(rightCircle);
        TwoDVector rightPointerVector = rightOrigin - leftOrigin;

        TwoDRightPointingLine line (leftOrigin, rightPointerVector);
        line.shortenByTwoCirclesRadiiAndExtraSpacing(
            leftCircle, rightCircle, extraSpacing
        );

        return line;
    }
};
