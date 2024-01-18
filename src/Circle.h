#pragma once

class Circle {
private:
    double xval, yval, rval;

public:
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
