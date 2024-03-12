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

};
