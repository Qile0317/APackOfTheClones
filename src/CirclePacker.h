#include <Rcpp.h>
#include <vector>
#include "../src/CircleNode.h"
#include "ProgressBar.h"

// general readability helpers

inline double sqr(double a) {
    return a * a;
}

inline double centre_dist(CircleNode& c) {
    return sqrt(sqr(c.x) + sqr(c.y));
}

inline std::pair<int, int> sentinel_pair() {
    return std::make_pair(-1, -1);
}

bool is_clear(std::pair<int, int>& p) {
    return (p.first == -1) && (p.second == -1);
}

class CirclePacker {
public:
    std::vector<CircleNode> data;
    int num_nodes;
    bool try_place;
    bool verbose;

    CirclePacker(
        const std::vector<double>& input_rad_vec,
        bool _try_place = true,
        bool _verbose = false
    ) {
        verbose = _verbose;
        start_progress_bar();

        try_place = _try_place;

        num_nodes = (int) input_rad_vec.size();
        data.resize(num_nodes);

        for (int i = 0; i < num_nodes; i++) {
            data[i] = CircleNode(input_rad_vec[i]);
        }
    }

    CirclePacker(
        const std::vector<CircleNode> final_circle_nodes,
        bool _try_place = true,
        bool _verbose = false
    ) {
        verbose = _verbose;
        start_progress_bar();

        try_place = _try_place;

        num_nodes = (int) final_circle_nodes.size();
        data = final_circle_nodes;
    }

    // main circle packing method, all helpers below.
    static Rcpp::List pack(
        std::vector<double> input_rad_vec,
        Rcpp::NumericVector centroid,
        double rad_decrease,
        bool _try_place,
        bool _verbose
    ) {
        CirclePacker packer = CirclePacker(input_rad_vec, _try_place, _verbose);

        if (packer.is_degenerate_case()) {
            return packer.handle_degenerate_cases(centroid, rad_decrease);
        }

        // for all nodes, fit node into cluster, check for overlap, and refit
        packer.place_starting_three();
        packer.init_boundary();

        int j = 3;
        while (j < packer.num_nodes) {
            int circ = packer.place_circle(j);
            packer.fit_circle(circ, j);
        }

        return packer.process_into_clusterlist(centroid, rad_decrease);
    }

    void progress_bar(int x, int max) {
        if (verbose) {
            ProgressBar::show(x, max);
        }
    }

    void start_progress_bar() {
        progress_bar(0, 1);
    }

    void finish_progress_bar() {
        progress_bar(1, 1);
    }

    // initialize node vector into the circular boundary linked list
    void init_boundary(int a = 0, int b = 2) {
        for (int i = a; i < b; i++) {
            data[i].nxt = i + 1;
            data[i + 1].prv = i;
        }
        data[b].nxt = a;
        data[b].prv = b - 1;
        data[a].prv = b;
    }

    // insert a node cn in-between c1 and c2 where inputs are indicies in data
    void insert_circle(int c1, int c2, int cn) {
        if (data[c1].nxt != c2 || data[c2].prv != c1) {
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
            curr = data[curr].nxt;
        }
        return dist;
    }

    // removes segment between c1 and c2
    void fwd_remove(int c1, int c2) {
        if (c1 == c2) { Rcpp::stop("Circles are the same."); }
        if (data[c1].nxt == c2) { Rcpp::stop("Circles are consecutive."); }

        data[c1].nxt = c2;
        data[c2].prv = c1;
    }

    // fit c3 tangent to c1 and c2, modifies data[c3].
    int fit_tang_circle(int c1, int c2, int c3) {
        double x1 = data[c1].x, x2 = data[c2].x;
        double y1 = data[c1].y, y2 = data[c2].y;
        double r1 = data[c1].rad, r2 = data[c2].rad;

        double r = data[c3].rad;
        double distance = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
        double invdist = 1 / distance;

        if (distance > (r1 + r2 + r + r)) {
            Rcpp::stop("Gap too large.");
        }

        double cos_sig = (x2 - x1) * invdist;
        double sin_sig = (y2 - y1) * invdist;
        double cos_gam = (sqr(distance)+sqr(r+r1)-sqr(r+r2))*0.5*invdist/(r+r1);
        double sin_gam = sqrt(1 - sqr(cos_gam));

        data[c3].x = x1 + ((r + r1) * (cos_sig * cos_gam - sin_sig * sin_gam));
        data[c3].y = y1 + ((r + r1) * (cos_sig * sin_gam + sin_sig * cos_gam));
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
        data[c1].x = -1 * data[c1].rad;
        data[c2].x = data[c2].rad;
        fit_tang_circle(c2, c1, c3); // the order is incredibly important

        double centroid_x = (data[c1].x + data[c2].x + data[c3].x) / 3;
        double centroid_y = (data[c1].y + data[c2].y + data[c3].y) / 3;

        data[c1].x -= centroid_x; data[c1].y -= centroid_y;
        data[c2].x -= centroid_x; data[c2].y -= centroid_y;
        data[c3].x -= centroid_x; data[c3].y -= centroid_y;
    }

    // finds the closest circle to the origin in the linked list containing c.
    int closest(int c) {
        int closest_c = c;
        int circ = data[c].nxt;
        while (circ != c) {
            if (centre_dist(data[closest_c]) > centre_dist(data[circ])) {
                closest_c = circ;
            }
            circ = data[circ].nxt;
        }
        return closest_c;
    }

    /*
    Locating the pair of successive circles, c, c.nxt, w/the following property:
    amongst all pairs of successive circles on the boundary, this pair minimizes
    distance from the center of d to the origin, when d is fitted tangent to
    this pair. Returns the index of the closest node.

    The function will modify some nodes in tang_circle_dist but when used later
    in circle_layout, the changes will be overwritten immediately after.
    */
    int closest_place(int c, int d) {
        int closest = c;
        int circ = data[c].nxt;
        while (circ != c) {
            double dist_closest = tang_circle_dist(
                closest, data[closest].nxt, d
            );
            double dist_circ = tang_circle_dist(circ, data[circ].nxt, d);

            if (dist_closest > dist_circ) {
                closest = circ;
            }
            circ = data[circ].nxt;
        }
        return closest;
    }

    // convenience function
    int place_circle(int j) {
        return (try_place) ? closest_place(j-1, j) : closest(j-1);
    }

    // check if two circle nodes overlap geometrically
    bool do_intersect(int c1, int c2) {
        double centre_distance = sqrt(
            sqr(data[c1].x - data[c2].x) + sqr(data[c1].y - data[c2].y)
        );
        double rad_sum = data[c1].rad + data[c2].rad;
        return centre_distance < rad_sum;
    }

    // convenience function for overlap_check
    int geod_dist(int Cm, int Cn, int C) {
        return std::min(fwd_dist(Cn, C), fwd_dist(C, Cm));
    }

    // construct a list of intersectors to node C from C_en to C_em exclusive
    std::vector<int> construct_obstruct_list(int C_en, int C_em, int C) {
        std::vector<int> obstruct;
        int circ = data[C_en].nxt;
        while (circ != C_em) {
            if (do_intersect(circ, C)) {
                obstruct.push_back(circ);
            }
            circ = data[circ].nxt;
        }
        return obstruct;
    }

    // overlap check of three circles
    std::pair<int, int> overlap_check(int Cm, int Cn, int C) {
        int C_em = Cm;
        int C_en = Cn;

        // collect circles that C intersects, and store indices in obstruct
        std::vector<int> obstruct = construct_obstruct_list(C_en, C_em, C);
        int n = obstruct.size();

        if (n == 0) {
            return sentinel_pair();
        }

        int nearest = obstruct[0];
        for (int i = 0; i < n; i++) {
            if (geod_dist(Cm, Cn, obstruct[i]) < geod_dist(Cm, Cn, nearest)) {
                nearest = obstruct[i];
            }
        }

        // if the distance is realized fwd, change C_en
        if (fwd_dist(Cn, nearest) <= fwd_dist(nearest, Cm)) {
            C_en = nearest;

        } else { // if distance is realized bkwd and not fwd, change C_em
            C_em = nearest;
        }

        return (C_em == Cm) && (C_en == Cn) ?
            sentinel_pair() : std::make_pair(C_em, C_en);
    }

    // helper for cases where there is only 1 or 2 input radii

    bool is_degenerate_case() {
        return ((num_nodes == 1) || (num_nodes == 2));
    }

    Rcpp::List handle_degenerate_cases(
        Rcpp::NumericVector& centroid,
        double rad_decrease
    ) {
        Rcpp::NumericVector x, y, r;
        double clrad;

        if (num_nodes == 1) {
            x = Rcpp::NumericVector::create(centroid[0]);
            y = Rcpp::NumericVector::create(centroid[1]);
            r = Rcpp::NumericVector::create(data[0].rad - rad_decrease);
            clrad = data[0].rad;

        } else if (num_nodes == 2) {
            x = Rcpp::NumericVector::create(
              (-1 * data[0].rad) + centroid[0], data[1].rad + centroid[0]
            );
            y = Rcpp::NumericVector::create(centroid[1], centroid[1]);
            r = Rcpp::NumericVector::create(
                data[0].rad - rad_decrease, data[1].rad - rad_decrease
            );
            clrad = 0.5 * (data[0].rad + data[1].rad);
        }

        Rcpp::List out = Rcpp::List::create(
            Rcpp::Named("x") = x,
            Rcpp::Named("y") = y,
            Rcpp::Named("rad") = r,
            Rcpp::Named("centroid") = centroid,
            Rcpp::Named("clRad") = clrad
        );

        finish_progress_bar();
        return out;
    }

    // returns the R list from a packed vector of circles for at least 3 circles
    Rcpp::List process_into_clusterlist(
        Rcpp::NumericVector& centroid,
        double rad_decrease
    ) {
        Rcpp::NumericVector x (num_nodes), y (num_nodes), rad (num_nodes);
        int max_x_ind = 0;

        for (int i = 0; i < num_nodes; i++) {
            x[i] = data[i].x + centroid[0];
            y[i] = data[i].y + centroid[1];
            rad[i] = data[i].rad - rad_decrease;

            if (x[i] > x[max_x_ind]) {
                max_x_ind = i;
            }
        }

        Rcpp::List out = Rcpp::List::create(
            Rcpp::Named("x") = x,
            Rcpp::Named("y") = y,
            Rcpp::Named("rad") = rad,
            Rcpp::Named("centroid") = centroid,
            Rcpp::Named("clRad") = x[max_x_ind] + rad[max_x_ind] - centroid[0]
        );

        finish_progress_bar();
        return out;
    }

    // fit a new circle (index j) adjacent to curr_circ and its neighbor
    void fit_circle(int curr_circ, int& j) {

        fit_tang_circle(curr_circ, data[curr_circ].nxt, j);
        std::pair<int, int> check = overlap_check(
            curr_circ, data[curr_circ].nxt, j
        );
        int cm = curr_circ, cn = data[curr_circ].nxt;

        while (!is_clear(check)) {
            cm = check.first;
            cn = check.second;

            fwd_remove(cm, cn);
            fit_tang_circle(cm, cn, j);
            check = overlap_check(cm, cn, j);
        }

        insert_circle(cm, cn, j);
        j++;

        if (j <= num_nodes) {
            progress_bar(j, num_nodes);
        }

    }

};
