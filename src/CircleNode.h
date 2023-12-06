class CircleNode {
public:
    double x;
    double y;
    double rad;
    int prv;
    int nxt;

    CircleNode() {
        x = 0;
        y = 0;
        rad = 0;
        prv = -1;
        nxt = -1;
    }

    CircleNode(double rad_val) {
        x = 0;
        y = 0;
        rad = rad_val;
        prv = -1;
        nxt = -1;
    }

    CircleNode(double x_val, double y_val, double rad_val) {
        x = x_val;
        y = y_val;
        rad = rad_val;
        prv = -1;
        nxt = -1;
    }

    bool operator==(CircleNode& other) {
        return (x == other.x)
        && (y == other.y)
        && (rad == other.rad)
        && (prv == other.prv)
        && (nxt == other.nxt);
    }
};
