#pragma once

class Circle {
private:
    double xval, yval, rval;

public:
    Circle() : xval(0), yval(0), rval(0) {}
    Circle(double x, double y, double r) : xval(x), yval(y), rval(r) {}
   
    bool hasOriginMoreLeftThan(Circle& other) {
        return x() < other.x();
    }

    // member var mutators

    Circle& increaseRad(double x) {
        return setRad(rval + x);
    }

    Circle& decreaseRad(double x) {
        return setRad(rval - x);
    }

    Circle& scaleRad(double scaleFactor) {
        return setRad(rval * scaleFactor);
    }

    // getters

    double x() const {
        return xval;
    }

    double y() const {
        return yval;
    }

    double rad() const {
        return rval;
    }

    // setters

    Circle& setX(double nx) {
        xval = nx;
        return *this;
    }

    Circle& setY(double ny) {
        yval = ny;
        return *this;
    }

    Circle& setRad(double nr) {
        rval = nr;
        return *this;
    }

};
