#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "data_structures.h"

using namespace std;

// initialize 3 nodes(ptrs) into the circular boundary linked list
void init_boundary(Node& c1, Node& c2, Node& c3) {
    c1.nxt = &c2; c1.prv = &c3;
    c2.nxt = &c3; c2.prv = &c1;
    c3.nxt = &c1; c3.prv = &c2;
}

// insert a node cn in-between c1 and c2
void insert_circle(Node& c1, Node& c2, Node& cn) {
    if (c1.nxt != &c2 || c2.prv != &c1) {
        throw invalid_argument("Two circles not adjacent");
    } else {
        c1.nxt = &cn;
        cn.prv = &c1;
        c2.prv = &cn;
        cn.nxt = &c2;
    }
}

// function to find "distance" between 2 elements of linked list
int fwd_dist(Node& c1, Node& c2) {
    int count = 0;
    Node* circ = &c1;
    while (circ != &c2) {
        count++;
        circ = circ->nxt;
    }
    return count;
}

// removes all nodes between c1 and c2 inlusive, and links everything else. Purposefully doesn't de-allocate the segment.
// this one simplies the original though doesnt unlink the removed segment. if there are bugs then this may be the cause.
void fwd_remove(Node& c1, Node& c2) {
    if (&c1 == &c2) { throw invalid_argument("Circles are the same."); }
    if (c1.nxt == &c2) { throw invalid_argument("Circles are consecutive."); }
    c1.prv->nxt = c2.nxt;
}

double centre_dist(Node& c) {
    return sqrt(c.x*c.x + c.y*c.y);
}

inline double sqr(double n) {
    return n*n;
}

// fit C3 tangent to C1 and C2, modified C3. Unsure about the return statement but works in R. probably will get errors later 
Node& fit_tang_circle(Node& C1, Node& C2, Node& C3) {
    double x1 = C1.x, x2 = C2.x;
    double y1 = C1.y, y2 = C2.y;
    double r1 = C1.rad, r2 = C2.rad;
  
    double r = C3.rad;
    double distance = sqrt(sqr(x1 - x2) + sqr(y1 - y2));

    if (distance > (r1 + r2 + r + r)) {
        throw invalid_argument("Gap too large.");
    }

    double cos_sig = (x2 - x1) / distance;
    double sin_sig = (y2 - y1) / distance;
    double cos_gam = (sqr(distance) + sqr(r + r1) - sqr(r + r2)) / (2 * distance * (r + r1));
    double sin_gam = sqrt(1 - sqr(cos_gam));

    C3.x = x1 + (r + r1) * (cos_sig * cos_gam - sin_sig * sin_gam);
    C3.y = y1 + (r + r1) * (cos_sig * sin_gam + sin_sig * cos_gam);
    return C3;
}

// place three mutually tangent circles to 0,0
void place_starting_three(Node& C1, Node& C2, Node& C3) {
    C1.x = -1 * C1.rad;
    C2.x = C2.rad;
    fit_tang_circle(C2, C1, C3); // the order is incredibly important. Also, surprised that there were no compiler/runtime errors though fit_tang_circle returns a ndoe reference.

    double centroid_x = (C1.x + C2.x + C3.x) / 3;
    double centroid_y = (C1.y + C2.y + C3.y) / 3;

    C1.x -= centroid_x;
    C2.x -= centroid_x;
    C3.x -= centroid_x;

    C1.y -= centroid_y;
    C2.y -= centroid_y;
    C3.y -= centroid_y;
}

// finds the closest circle to the origin in the linked list containing c. needs testing
Node& closest(Node& c) {
    Node& closest = c;
    Node& circ = *(c.nxt);

    while (&circ != &c) {
        if (centre_dist(closest) > centre_dist(circ)) {
            closest = circ;
        }
        circ = *(circ.nxt);
    }
    return closest;
}

// unfinished