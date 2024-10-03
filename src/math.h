#pragma once

#include "Circle.h"

#include <cmath>

#define TAU 6.28318530717958647692

inline double sqr(double x) {
    return x * x;
}

inline double eucDist(double x1, double x2) {
    return std::sqrt(sqr(x1) + sqr(x2));
}

// TODO test these

inline double normalizeAngle(double radianAngle) {
    return std::abs(radianAngle) <= TAU ? radianAngle
        : fmod(radianAngle, TAU);
}

inline double normalizeAngleAtan2(double radianAngle) {
    radianAngle = normalizeAngle(radianAngle);
    if (std::abs(radianAngle) <= M_PI) return radianAngle;
    return std::abs(radianAngle) - ((1 + (radianAngle > 0)) * M_PI);
}

// representation of vector in R^2. getMagnitude is guaranteed to be
// always positive, and getDirection is guaranteed to be in (-pi, pi)
class TwoDVector {
protected:

    double magnitude;
    double direction;

public:

    // constructors

    TwoDVector() : magnitude(0), direction(0) {}
   
    TwoDVector(double mag, double dir) : magnitude(mag), direction(dir) {}
   
    TwoDVector(const TwoDVector& v) :
        magnitude(v.getMagnitude()), direction(v.getDirection()) {}

    // factory constructors

    static TwoDVector createFromXY(double x, double y) {
        return TwoDVector(eucDist(x, y), std::atan2(y, x));
    }

    static TwoDVector createCircleOrigin(const Circle& circle) {
        return createFromXY(circle.x(), circle.y());
    }

    // copy the values of another vector

    TwoDVector& copyFrom(const TwoDVector& other) {
        magnitude = other.getMagnitude();
        direction = other.getDirection();
        return *this;
    }

    // 2d vector maths

    TwoDVector operator+(const TwoDVector& other) const {
        return createFromXY(getX() + other.getX(), getY() + other.getY());
    }

    TwoDVector& operator+=(const TwoDVector& other) {
        TwoDVector result = *this + other;
        return copyFrom(result);
    }

    TwoDVector operator-() const {
        return createFromXY(-getX(), -getY()); // can also use reverse direction
    }

    TwoDVector operator-(const TwoDVector& other) const {
        return *this + (-other);
    }

    bool operator==(const TwoDVector& other) const {
        return getX() == other.getX() && getY() == other.getY();
    }

    // component coordinate mutation math

    TwoDVector& increaseXYComponent(double xval, double yval) {
        return decreaseXYComponent(-xval, -yval);
    }

    TwoDVector& decreaseXYComponent(double xval, double yval) {
        return decreaseXComponent(xval).decreaseYComponent(yval);
    }

    TwoDVector& decreaseXComponent(double val) {
        return copyFrom(createFromXY(getX() - val, getY()));
    }

    TwoDVector& decreaseYComponent(double val) {
        return copyFrom(createFromXY(getX(), getY() - val));
    }

    // polar coordinate mutation math

    TwoDVector& increaseMagnitude(double val) {
        return setMagnitude(getMagnitude() + val);
    }

    TwoDVector& decreaseMagnitude(double val) {
        return increaseMagnitude(-val);
    }

    TwoDVector& scaleMagnitude(double val) {
        return setMagnitude(getMagnitude() * val);
    }

    TwoDVector& reverseDirection() {
        return setDirection(direction + M_PI);
    }

    // setters

    TwoDVector& setMagnitude(double val) {
        magnitude = std::abs(val);
        return val > 0 ? *this : reverseDirection();
    }

    TwoDVector& setDirection(double val) {
        direction = normalizeAngleAtan2(val);
        return *this;
    }

    // getters

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

};

class TwoDLine {
protected:

    TwoDVector leftPoint;
    TwoDVector rightPoint;

public:

    TwoDLine(const TwoDVector& left, const TwoDVector& right) {
        leftPoint = TwoDVector(left);
        rightPoint = TwoDVector(right);
    }

    TwoDLine& shortenByTwoCirclesRadii(double leftRadius, double rightRadius) {
        return shortenLeftByCircleRadius(leftRadius)
            .shortenRightByCircleRadius(rightRadius);
    }

    TwoDLine& shortenLeftByCircleRadius(double val) {
        leftPoint.copyFrom(
            rightPoint + (leftPoint - rightPoint).decreaseMagnitude(val)
        );
        return *this;
    }

    TwoDLine& shortenRightByCircleRadius(double val) {
        rightPoint.copyFrom(
            leftPoint + (rightPoint - leftPoint).decreaseMagnitude(val)
        );
        return *this;
    }

    int matchLeftCircleIndex(std::vector<Circle>& circles, int i, int j) {
        return leftPoint == TwoDVector::createCircleOrigin(circles[i]) ? i : j;
    }

    int matchRightCircleIndex(std::vector<Circle>& circles, int i, int j) {
        return i + j - matchLeftCircleIndex(circles, i, j);
    }

    // getters

    TwoDVector getLeftPoint() const {
        return leftPoint;
    }

    TwoDVector getRightPoint() const {
        return rightPoint;
    }

    double getLeftX() const {
        return getLeftPoint().getX();
    }

    double getLeftY() const {
        return getLeftPoint().getY();
    }

    double getRightX() const {
        return getRightPoint().getX();
    }

    double getRightY() const {
        return getRightPoint().getY();
    }

};

class TwoDLineFactory {
public:

    static TwoDLine createCircleLinkingLineWithSpacing(
        Circle& c1, Circle& c2, double extraSpacing
    ) {

        Circle& leftCircle = c1, rightCircle = c2;
        if (c2.hasOriginMoreLeftThan(c1)) std::swap(leftCircle, rightCircle);

        TwoDLine circleLinkingLine = TwoDLine(
            TwoDVector::createCircleOrigin(leftCircle),
            TwoDVector::createCircleOrigin(rightCircle)
        );

        return circleLinkingLine.shortenByTwoCirclesRadii(
            leftCircle.rad() + extraSpacing, rightCircle.rad() + extraSpacing
        );
    }

};
