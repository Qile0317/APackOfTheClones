/*
#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
*/
 
// node of circular linked list where each node is a circle
class Node {
public:
  double x;      
  double y;      
  double rad;    
  Node* prv; 
  Node* nxt; 
  
  Node(double x_val, double y_val, double rad_val) {
    x = x_val;
    y = y_val;
    rad = rad_val;
    prv = nullptr;
    nxt = nullptr;
  }
  
  // node printer for debugging
  friend std::ostream& operator<<(std::ostream& os, const Node& node) {
    os << "x: " << node.x << ", y: " << node.y << ", rad: " << node.rad;
    return os;
  }
  
  // traverse and print method for debugging
  void traverse() {
    Node* curr = this;
    while (curr != nullptr && curr->nxt != this) {
      std::cout << *curr << std::endl;
      curr = curr->nxt;
    }
    std::cout << *curr << std::endl;
  }
};

// initialize node vector into the circular boundary linked list
void init_boundary(std::vector<Node*> nodes) {
  int n = nodes.size() - 1;
  for (int i = 0; i < n; i++) {
    nodes[i+1]->prv = nodes[i];
    nodes[i]->nxt = nodes[i+1];
  }
  nodes[n]->nxt = nodes[0];
  nodes[0]->prv = nodes[n];
}

// insert a node cn in-between c1 and c2
void insert_circle(Node& c1, Node& c2, Node& cn) {
  if (c1.nxt != &c2 || c2.prv != &c1) {
    Rcpp::stop("Two circles not adjacent");
  } else {
    c1.nxt = &cn;
    cn.prv = &c1;
    c2.prv = &cn;
    cn.nxt = &c2;
  }
}

// function to find forward "distance" between 2 elements of linked list
int fwd_dist(Node& c1, Node& c2) {
  int count = 0;
  Node* circ = &c1;
  while (circ != &c2) {
    count++;
    circ = circ->nxt;
  }
  return count;
}

// removes segment between c1 and c2 as one moves forwards
void fwd_remove(Node& c1, Node& c2) {
  if (&c1 == &c2) { Rcpp::stop("Circles are the same."); }
  if (c1.nxt == &c2) { Rcpp::stop("Circles are consecutive."); }
  
  c1.nxt = &c2;
  c2.prv= &c1;
  
  /*
  Node* circ = c1.nxt;
  while (circ != &c2) {
    circ->prv->nxt = circ->nxt;
    circ->nxt->prv = circ->prv;
    circ = circ->nxt;
  }
  */
}

inline double centre_dist(Node& c) {
  return sqrt(sqr(c.x) + sqr(c.y));
}

// fit C3 tangent to C1 and C2, modified C3. Unsure about the return statement but works in R. probably will get errors later 
Node& fit_tang_circle(Node& C1, Node& C2, Node& C3) {
  double x1 = C1.x, x2 = C2.x;
  double y1 = C1.y, y2 = C2.y;
  double r1 = C1.rad, r2 = C2.rad;
  
  double r = C3.rad;
  double distance = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
  
  if (distance > (r1 + r2 + r + r)) {
    Rcpp::stop("Gap too large.");
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
// fits C3 such that C1,C2,C3 are arranged counterclockwise
void place_starting_three(Node& C1, Node& C2, Node& C3) {
  C1.x = -1 * C1.rad;
  C2.x = C2.rad;
  fit_tang_circle(C2, C1, C3); // the order is incredibly important. Also, surprised that there were no compiler/runtime errors though fit_tang_circle returns a ndoe reference.
  
  double centroid_x = (C1.x + C2.x + C3.x) / 3;
  double centroid_y = (C1.y + C2.y + C3.y) / 3;
  
  C1.x -= centroid_x; C1.y -= centroid_y;
  C2.x -= centroid_x; C2.y -= centroid_y;
  C3.x -= centroid_x; C3.y -= centroid_y;
}

// finds the closest circle to the origin in the linked list containing c. Returns a pointer to the node.
Node* closest(Node& c) {
  Node* closest = &c;
  Node* circ = c.nxt;
  while (circ != &c) {
    if (centre_dist(*(closest)) > centre_dist(*(circ))) {
      closest = circ;
    }
    circ = circ->nxt;
  }
  return closest;
}

/* 
 Locating the pair of successive circles, c, c.nxt, with the following property: 
 amongst all pairs of successive circles on the boundary, this pair minimizes distance from the center of d to the origin, when d is
 fitted tangent to this pair. Returns a pointer to the output node.
 */
Node* closest_place(Node& c1, Node& c2) {
  Node* closest = &c1;
  Node* circ = c1.nxt;
  while (circ != &c1) {
    if (centre_dist(fit_tang_circle(*(closest), *(closest->nxt), c2)) > centre_dist(fit_tang_circle(*(circ), *(circ->nxt), c2))) {
      closest = circ;
    }
    circ = circ->nxt;
  }
  return closest;
}

// check if two circle nodes overlap geometrically
inline bool do_intersect(Node& c1, Node& c2) {
  return sqrt(sqr(c1.x - c2.x) + sqr(c1.y - c2.y)) < ((c1.rad + c2.rad));
}

inline int geod_dist(Node& Cm, Node& Cn, Node& C) {
  return std::min(fwd_dist(Cn, C), fwd_dist(C, Cm));
}

// overlap check of three circles and returns either a vector of two pointers to nodes or nothing
// unfinished
std::pair<Node*, Node*> overlap_check(Node& Cm, Node& Cn, Node& C) {
  Node* C_em = &Cm;
  Node* C_en = &Cn;
  std::vector<Node*> obstruct;
  
  // collect circles that C intersects (if any) by adding intersectors to 'obstruct'
  Node* circ = Cn.nxt;
  while (circ != &Cm) {
    if (do_intersect(*(circ), C)) { // not sure if doing *(circ) here still makes reference obj
      obstruct.push_back(circ);
    }
    circ = circ->nxt;
  }
  
  // find the one closest to {Cm, Cn}, (distance is in number of steps)
  int LenObs = obstruct.size(); 
  if (LenObs > 0) { 
    Node* nearest = obstruct[0];
    for (int i = 0; i < LenObs; i++) {
      if (geod_dist(Cm, Cn, *(obstruct[i])) < geod_dist(Cm, Cn, *(nearest))) {
        nearest = obstruct[i];
      }
    }
    
    if (fwd_dist(Cn, *(nearest)) <= fwd_dist(*(nearest), Cm)) { // if the distance is realized fwd, change C_en
      C_en = nearest;
    } else { // if distance is realized bkwd and not fwd, change C_em
      C_em = nearest;
    }
  }
  if ((C_em == &Cm) && (C_en == &Cn)) {
    return std::make_pair(nullptr, nullptr);
  }
  return std::make_pair(C_em, C_en);
}

inline bool is_degenerate_case(int lenCirc) {
  bool res = (lenCirc == 1) || (lenCirc == 2);
  return res;
}

// unfinished
Rcpp::List handle_degenerate_cases(
  int lenCirc,
  std::vector<Node>& circles,
  std::pair<double, double>& centroid,
  double rad_scale_factor,
  bool verbose
) {
  if (lenCirc == 1) {
    return Rcpp::List::create();
  }
  return Rcpp::List::create();
}

// [[Rcpp::export]]
double estimate_rad(
    std::vector<double> x_vals,
    std::vector<double> rad_vals,
    double centroid_x
) {
  double max_x = 0;
  int max_ind = 0, n = x_vals.size();
  for (int i = 0; i < n; i++) {
    if (x_vals[i] > max_x) {
      max_x = x_vals[i];
      max_ind += i - max_ind;
    }
  }
  double result_num = max_x + rad_vals[max_ind] - centroid_x;
  return result_num;
}