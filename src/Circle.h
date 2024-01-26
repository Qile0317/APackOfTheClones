#pragma once

class Circle {
private:
    double xval, yval, rval;

public:
    Circle() : xval(0), yval(0), rval(0) {}
    Circle(double x, double y, double r) : xval(x), yval(y), rval(r) {}
    
    double x() {
        return xval;
    }

    double y() {
        return yval;
    }

    double rad() {
        return rval;
    }

    bool hasOriginMoreLeftThan(Circle& other) {
        return xval < other.xval;
    }
};
