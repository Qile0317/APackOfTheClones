#pragma once

class Circle {
private:
    double xval, yval, rval;

public:
    Circle() : xval(0), yval(0), rval(0) {}
    Circle(double x, double y, double r) : xval(x), yval(y), rval(r) {}
    
    int x() {
        return xval;
    }

    int y() {
        return yval;
    }

    int rad() {
        return rval;
    }
};
