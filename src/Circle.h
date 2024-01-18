class Circle {
private:
    double x, y, r;

public:
    Circle(double x, double y, double r) : x(x), y(y), r(r) {}
    
    int x() {
        return x;
    }

    int y() {
        return y;
    }

    int rad() {
        return r;
    }
};
