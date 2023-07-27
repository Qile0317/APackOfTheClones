// script to pack a vector of radii into a circular cluster of circles

class Node {
public:
    double x;
    double y;
    double rad;
    int prv;
    int nxt;

    Node(double x_val = 0, double y_val = 0, double rad_val = 0) {
        x = x_val;
        y = y_val;
        rad = rad_val;
        prv = -1;
        nxt = -1;
    }

    bool operator==(Node& other) {
        return (x == other.x)
        && (y == other.y)
        && (rad == other.rad)
        && (prv == other.prv)
        && (nxt == other.nxt);
    }
};

// functions for readability
inline double centre_dist(Node& c) {
    return sqrt(sqr(c.x) + sqr(c.y));
}

inline bool is_clear(std::pair<int, int>& p) {
    return (p.first == -1) && (p.second == -1);
}

// macros for slightly better syntax for accessing Node attributes in member functions
#define _x(a) (data[a].x)
#define _y(a) (data[a].y)
#define _rad(a) (data[a].rad)
#define nxt(a) (data[a].nxt)
#define prv(a) (data[a].prv)

class NodeVector {
public:
    std::vector<Node> data;
    int num_nodes;

    // Constructors
    NodeVector(const std::vector<Node> nodes) {
        num_nodes = nodes.size();
        data = nodes;
    }

    NodeVector(const std::vector<double>& input_rad_vec, int _num_nodes = 0) {
        if (_num_nodes == 0) {
            num_nodes = (int)input_rad_vec.size();
        }
        data.resize(num_nodes);

        for (int i = 0; i < num_nodes; i++) {
            data[i] = Node(0, 0, input_rad_vec[i]);
        }
    }

    // functions for cpp_circle_layout (exported function at end of script)

    // initialize node vector into the circular boundary linked list
    void init_boundary(int a = 0, int b = 2) {
        for (int i = a; i < b; i++) {
            data[i].nxt = i+1;
            data[i+1].prv = i;
        }
        data[b].nxt = a;
        data[b].prv = b-1;
        data[a].prv = b;
    }

    // insert a node cn in-between c1 and c2 where inputs are indicies in data
    void insert_circle(int c1, int c2, int cn) {
        if (nxt(c1) != c2 || prv(c2) != c1) {
            Rcpp::stop("Two circles not adjacent");
            return;
        }
        data[c1].nxt = cn;
        data[cn].prv = c1;
        data[c2].prv = cn;
        data[cn].nxt = c2;
    }

    // function to find forward "distance" between 2 elements of linked list
    int fwd_dist(int c1, int c2) {
        int curr = c1;
        int dist = 0;
        while (curr != c2) {
            dist++;
            curr = nxt(curr);
        }
        return dist;
    }

    // removes segment between c1 and c2, Rcpp::Stop technically not needed as users wont access this
    void fwd_remove(int c1, int c2) {
        if (c1 == c2) { Rcpp::stop("Circles are the same."); return; }
        if (nxt(c1) == c2) { Rcpp::stop("Circles are consecutive."); return; }

        data[c1].nxt = c2;
        data[c2].prv = c1;
    }

    // fit c3 tangent to c1 and c2, modifies data[c3].
    int fit_tang_circle(int c1, int c2, int c3) {
        double x1 = _x(c1), x2 = _x(c2);
        double y1 = _y(c1), y2 = _y(c2);
        double r1 = _rad(c1), r2 = _rad(c2);

        double r = _rad(c3);
        double distance = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
        double invdist = 1/distance;

        if (distance > (r1 + r2 + r + r)) { // probably not needed
            Rcpp::stop("Gap too large.");
        }

        double cos_sig = (x2 - x1) * invdist;
        double sin_sig = (y2 - y1) * invdist;
        double cos_gam = (sqr(distance)+sqr(r+r1) - sqr(r+r2)) * 0.5*invdist / (r+r1);
        double sin_gam = sqrt(1 - sqr(cos_gam));

        _x(c3) = x1 + ((r + r1) * (cos_sig * cos_gam - sin_sig * sin_gam));
        _y(c3) = y1 + ((r + r1) * (cos_sig * sin_gam + sin_sig * cos_gam));
        return c3;
    }

    // convenience function for closest_place
    // It modifies the objects but when its used later, it "cancels out"
    double tang_circle_dist(int c1, int c2, int c3) {
        return centre_dist(data[fit_tang_circle(c1, c2, c3)]);
    }

    // place three mutually tangent circles to 0,0
    // fits c3 such that c1,c2,c3 are arranged counterclockwise
    void place_starting_three(int c1 = 0, int c2 = 1, int c3 = 2) {
        data[c1].x = -1 * _rad(c1);
        data[c2].x = _rad(c2);
        fit_tang_circle(c2, c1, c3); // the order is incredibly important

        double centroid_x = (_x(c1) + _x(c2) + _x(c3)) / 3;
        double centroid_y = (_y(c1) + _y(c2) + _y(c3)) / 3;

        data[c1].x -= centroid_x; data[c1].y -= centroid_y;
        data[c2].x -= centroid_x; data[c2].y -= centroid_y;
        data[c3].x -= centroid_x; data[c3].y -= centroid_y;
    }

    // finds the closest circle to the origin in the linked list containing c.
    int closest(int c) {
        int closest_c = c;
        int circ = nxt(c);
        while (circ != c) {
            if (centre_dist(data[closest_c]) > centre_dist(data[circ])) {
                closest_c = circ;
            }
            circ = nxt(circ);
        }
        return closest_c;
    }

    /*
    Locating the pair of successive circles, c, c.nxt, with the following property:
    amongst all pairs of successive circles on the boundary, this pair minimizes
    distance from the center of d to the origin, when d is fitted tangent to this
    pair. Returns the index of the closest node

    The function will modify some nodes in tang_circle_dist but when used later in
    circle_layout, the changes will be overwritten immediately after
    */
    int closest_place(int c1, int c2) {
        int closest = c1;
        int circ = nxt(c1);
        while (circ != c1) {
            double dist_closest = tang_circle_dist(closest, data[closest].nxt, c2);
            double dist_circ = tang_circle_dist(circ, nxt(circ), c2);

            if (dist_closest > dist_circ) {
                closest = circ;
            }
            circ = nxt(circ);
        }
        return closest;
    }

    // convinience function
    int place_circle(int j, bool try_place) {
        return (try_place) ? closest_place(j-1, j) : closest(j-1);
    }

    // check if two circle nodes overlap geometrically
    bool do_intersect(int c1, int c2) {
        double centre_distance = sqrt(sqr(_x(c1) - _x(c2)) + sqr(_y(c1) - _y(c2)));
        double rad_sum = _rad(c1) + _rad(c2);
        return centre_distance < rad_sum;
    }

    // convenience function for overlap_check
    int geod_dist(int Cm, int Cn, int C) {
        return std::min(fwd_dist(Cn, C), fwd_dist(C, Cm));
    }

    // construct a list of intersectors to node C from C_en to C_em exclusive
    std::vector<int> construct_obstruct_list(int C_en, int C_em, int C) {
        std::vector<int> obstruct;
        int circ = nxt(C_en);
        while (circ != C_em) {
            if (do_intersect(circ, C)) {
                obstruct.push_back(circ);
            }
            circ = nxt(circ);
        }
        return obstruct;
    }

    // overlap check of three circles - seems to be O(length(nodes)^2)
    std::pair<int, int> overlap_check(int Cm, int Cn, int C) {
        int C_em = Cm;
        int C_en = Cn;

        // collect circles that C intersects, if any, and store indicies in obstruct
        std::vector<int> obstruct = construct_obstruct_list(C_en, C_em, C);
        int n = obstruct.size();

        if (n == 0) {
            return std::make_pair(-1, -1);
        }

        int nearest = obstruct[0];
        for (int i = 0; i < n; i++) {
            if (geod_dist(Cm, Cn, obstruct[i]) < geod_dist(Cm, Cn, nearest)) {
                nearest = obstruct[i];
            }
        }

        if (fwd_dist(Cn, nearest) <= fwd_dist(nearest, Cm)) { // if the distance is realized fwd, change C_en
            C_en = nearest;
        } else { // if distance is realized bkwd and not fwd, change C_em
            C_em = nearest;
        }


        if ((C_em == Cm) && (C_en == Cn)) {
            return std::make_pair(-1, -1);
        }
        return std::make_pair(C_em, C_en);
    }

    bool is_degenerate_case() {
        bool res = (num_nodes == 1) || (num_nodes == 2);
        return res;
    }

    Rcpp::List handle_degenerate_cases(
        Rcpp::NumericVector& centroid,
        double rad_decrease,
        bool verbose = false
    ) {
        Rcpp::NumericVector x, y, r;
        double clrad;

        if (num_nodes == 1) {
            x = Rcpp::NumericVector::create(centroid[0]);
            y = Rcpp::NumericVector::create(centroid[1]);
            r = Rcpp::NumericVector::create(_rad(0) - rad_decrease);
            clrad = _rad(0);

        } else if (num_nodes == 2) {
            x = Rcpp::NumericVector::create(
              (-1 * _rad(0)) + centroid[0], _rad(1) + centroid[0]
            );
            y = Rcpp::NumericVector::create(centroid[1], centroid[1]);
            r = Rcpp::NumericVector::create(
                _rad(0) - rad_decrease, _rad(1) - rad_decrease
            );
            clrad = 0.5 * (_rad(0) + _rad(1));
        }

        if (verbose) {progress_bar(1, 1);}
        return Rcpp::List::create(
            _["x"] = x, _["y"] = y, _["rad"] = r,
            _["centroid"] = centroid, _["clRad"] = clrad
        );
    }

    // return the R list from a packed vector of circles for at least 3 circles
    Rcpp::List process_into_clusterlist(
        Rcpp::NumericVector& centroid,
        double rad_decrease,
        bool verbose = false
    ) {
        Rcpp::NumericVector x (num_nodes), y (num_nodes), rad (num_nodes);
        int max_x_index = 0;

        for (int i = 0; i < num_nodes; i++) {
            x[i] = _x(i) + centroid[0];
            y[i] = _y(i) + centroid[1];
            rad[i] = _rad(i) - rad_decrease;

            if (x[i] > x[max_x_index]) {
                max_x_index = i;
            }
        }

        if(verbose) {
            progress_bar(1, 1);
        }
        return Rcpp::List::create(
            _["x"] = x, _["y"] = y, _["rad"] = rad, _["centroid"] = centroid,
            _["clRad"] = x[max_x_index] + rad[max_x_index] - centroid[0]
        );
    }

    // fit a new circle (index j) adjacent to curr_circ and its neighbour
    int fit_circle(int curr_circ, int& j, bool verbose) {

        fit_tang_circle(curr_circ, nxt(curr_circ), j);
        std::pair<int, int> check = overlap_check(curr_circ, nxt(curr_circ), j);
        int cm = curr_circ, cn = nxt(curr_circ);

        while (!is_clear(check)) {
            cm = check.first;
            cn = check.second;

            fwd_remove(cm, cn);
            fit_tang_circle(cm, cn, j);
            check = overlap_check(cm, cn, j);
        }

        insert_circle(cm, cn, j);
        j++;

        if (verbose && (j <= num_nodes)) {
            progress_bar(j, num_nodes);
        }
        return j; // probably unessecary?
    }
};

// [[Rcpp::export]]
Rcpp::List cpp_circle_layout(
    std::vector<double> input_rad_vec,
    Rcpp::NumericVector centroid,
    double rad_decrease = 0,
    bool try_place = false,
    bool verbose = true
) {
    if(verbose) {progress_bar(0, 1);}

    NodeVector nodes = NodeVector(input_rad_vec);

    if (nodes.is_degenerate_case()) {
        return nodes.handle_degenerate_cases(centroid, rad_decrease, verbose);
    }

    // for all nodes, fit each node into cluster, check for overlap, and refit
    nodes.place_starting_three();
    nodes.init_boundary();

    int j = 3;
    while (j < nodes.num_nodes) {
        int circ = nodes.place_circle(j, try_place);
        j = nodes.fit_circle(circ, j, verbose);
    }

    return nodes.process_into_clusterlist(centroid, rad_decrease, verbose);
}
