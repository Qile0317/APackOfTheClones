// node of circular linked list where each node is a circle
class Node {
  public:
    double x;      
    double y;      
    double rad;    
    Node* prv; 
    Node* nxt; 
  
    Node(double x_val = 0, double y_val = 0, double rad_val = 0) {
      x = x_val;
      y = y_val;
      rad = rad_val;
      prv = nullptr;
      nxt = nullptr;
    }
  
    bool operator==(Node& other) {
      return (x == other.x)
        && (y == other.y)
        && (rad == other.rad)
        && (prv == other.prv)
        && (nxt == other.nxt);
    }
  
  private:
    friend std::ostream& operator<<(std::ostream& os, Node& node) {
      os << "x: " << node.x << ", y: " << node.y << ", rad: " << node.rad;
      return os;
    }
    
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
  Node* target_circ = &c2;
  while (circ != target_circ) {
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
  c2.prv = &c1;
  
  /*
  Node* circ = c1.nxt, target_circ = &c2;;
  while (circ != target_circ) {
    circ->prv->nxt = circ->nxt;
    circ->nxt->prv = circ->prv;
    circ = circ->nxt;
  }
  */
}

// euclidean distance of circle center to (0, 0)
double centre_dist(Node& c) {
  return sqrt(sqr(c.x) + sqr(c.y));
}

// fit C3 tangent to C1 and C2, modifies C3. 
Node& fit_tang_circle(Node& C1, Node& C2, Node& C3) {
  double x1 = C1.x, x2 = C2.x;
  double y1 = C1.y, y2 = C2.y;
  double r1 = C1.rad, r2 = C2.rad;
  
  double r = C3.rad;
  double distance = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
  double invdist = 1/distance; // could use the fast inverse square root here
  
  if (distance > (r1 + r2 + r + r)) {
    Rcpp::stop("Gap too large.");
  }
  
  double cos_sig = (x2 - x1) * invdist;
  double sin_sig = (y2 - y1) * invdist;
  double cos_gam = (sqr(distance)+sqr(r+r1) - sqr(r+r2)) * 0.5*invdist / (r+r1);
  double sin_gam = sqrt(1 - sqr(cos_gam));
  
  C3.x = x1 + (r + r1) * (cos_sig * cos_gam - sin_sig * sin_gam);
  C3.y = y1 + (r + r1) * (cos_sig * sin_gam + sin_sig * cos_gam);
  return C3;
}

// convenience function for closest_place
// It modifies the objects but when its used later, it "cancels out"
double tang_circle_dist(Node& C1, Node& C2, Node& C3) {
  return centre_dist(fit_tang_circle(C1, C2, C3));
}

// place three mutually tangent circles to 0,0
// fits C3 such that C1,C2,C3 are arranged counterclockwise
void place_starting_three(Node& C1, Node& C2, Node& C3) {
  C1.x = -1 * C1.rad;
  C2.x = C2.rad;
  fit_tang_circle(C2, C1, C3); // the order is incredibly important
  
  double centroid_x = (C1.x + C2.x + C3.x) / 3;
  double centroid_y = (C1.y + C2.y + C3.y) / 3;
  
  C1.x -= centroid_x; C1.y -= centroid_y;
  C2.x -= centroid_x; C2.y -= centroid_y;
  C3.x -= centroid_x; C3.y -= centroid_y;
}

// finds the closest circle to the origin in the linked list containing c. Returns a pointer to the node.
Node& closest(Node& c) {
  Node* closest = &c;
  Node* circ = c.nxt;
  while (circ != &c) {
    if (centre_dist(*(closest)) > centre_dist(*(circ))) {
      closest = circ;
    }
    circ = circ->nxt;
  }
  return *(closest);
}

/* 
 Locating the pair of successive circles, c, c.nxt, with the following property: 
 amongst all pairs of successive circles on the boundary, this pair minimizes
 distance from the center of d to the origin, when d is fitted tangent to this
 pair. Returns a pointer to the output node. 
 
 The function will modify some nodes but when used later in circle_layout, the
 changes will be overwritten immediately after
 */
Node& closest_place(Node& c1, Node& c2) {
  Node* closest = &c1;
  Node* circ = c1.nxt;
  while (circ != &c1) {
    double dist_closest = tang_circle_dist(*(closest),*(closest->nxt),c2);
    double dist_circ = tang_circle_dist(*(circ), *(circ->nxt), c2);
    
    if (dist_closest > dist_circ) {
      closest = circ;
    }
    circ = circ->nxt;
  }
  return *(closest);
}

// check if two circle nodes overlap geometrically
bool do_intersect(Node& c1, Node& c2) {
  return sqrt(sqr(c1.x - c2.x) + sqr(c1.y - c2.y)) < ((c1.rad + c2.rad));
}

// convenience function for overlap_check
int geod_dist(Node& Cm, Node& Cn, Node& C) {
  return std::min(fwd_dist(Cn, C), fwd_dist(C, Cm));
}

// helper function for overlap_check to create a sentinel node pointer pair
std::pair<Node*, Node*> make_sentinel_pair() {
  Node* nullPtr = nullptr;
  std::pair<Node*, Node*> SENTINEL_PAIR = std::make_pair(nullPtr, nullPtr);
  return SENTINEL_PAIR;
}

// helper function for overlap check to check if a pair is sentinel (clear)
// I actually think this probably takes advantage of a "bug" since to my
// understanding the pointers from make_sentinel_pair are null due to automatic
// memory deallocation but if it works it works :/
bool is_clear(std::pair<Node*, Node*>& inp) {
  return (inp.first == nullptr) && (inp.second == nullptr);
}

// overlap check of three circles and returns either a vector of two pointers
// or the sentinel node pair if check is "clear"
// Problem: assumes valid inputs where iterating over Cn will actually get to Cm
// will be infinite loop if not
std::pair<Node*, Node*> overlap_check(Node& Cm, Node& Cn, Node& C) {
  Node* C_em = &Cm;
  Node* C_en = &Cn;
  std::vector<Node*> obstruct;
  
  Rcpp::Rcout << "overlap check: initialized vars, starting while loop of death\n";
  // collect circles that C intersects (if any) by adding intersectors to 'obstruct'
  Node* circ = Cn.nxt;
  Rcpp::Rcout << "starting circ: " << circ << " C_em: " << C_em << "\n";
  while (circ != C_em) {
    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "-------------------------\ncirc: " << circ << "\ncirc.nxt: " << circ->nxt << "\n";

    Node& curr_circ = *(circ);
    Rcpp::Rcout << "intersect with C?: " << do_intersect(curr_circ, C) << "\n";
    if (do_intersect(curr_circ, C)) {
      obstruct.push_back(circ);
      Rcpp::Rcout << obstruct.size() << "\n";
    }
    circ = circ->nxt;
  }
  
  Rcpp::Rcout << "obstruct list construction complete\n";
  
  // find the one closest to {Cm, Cn}, (distance is in number of steps)
  int n = obstruct.size(); 
  if (n > 0) { 
    Node& nearest = *(obstruct[0]);
    for (int i = 0; i < n; i++) {
      Node& curr_obstruct_node = *(obstruct[i]);
      if (geod_dist(Cm, Cn, curr_obstruct_node) < geod_dist(Cm, Cn, nearest)) {
        nearest = curr_obstruct_node;
      }
    }
    
    if (fwd_dist(Cn, nearest) <= fwd_dist(nearest, Cm)) { // if the distance is realized fwd, change C_en
      C_en = &nearest;
    } else { // if distance is realized bkwd and not fwd, change C_em
      C_em = &nearest;
    }
  }
  Rcpp::Rcout << "pretty much done\n";
  if ((C_em == &Cm) && (C_en == &Cn)) {
    return make_sentinel_pair();
  }
  return std::make_pair(C_em, C_en);
}

// functions for handling "degenerate cases" when theres only 1 or 2 nodes in
// the radii list
bool is_degenerate_case(int n) {
  bool res = (n == 1) || (n == 2);
  return res;
}

// needs testing
Rcpp::List handle_degenerate_cases(
    int lenCirc,
    std::vector<Node>& circles, // not sure if node should be referenced
    Rcpp::NumericVector& centroid,
    double rad_scale_factor,
    bool verbose
) {
  Rcpp::NumericVector x, y, r;
  double clrad;
  
  if (lenCirc == 1) {
    x = {centroid[0]};
    y = {centroid[1]};
    r = {circles[0].rad * rad_scale_factor};
    clrad = r[0];
  } else {
    x = {(-1*circles[0].rad) + centroid[0], circles[1].rad + centroid[0]};
    y = {centroid[1], centroid[1]};
    r = {circles[0].rad * rad_scale_factor,
         circles[1].rad * rad_scale_factor
        };
    clrad = 0.5 * (r[0] + r[1]);
  }
  if (verbose) {progress_bar(1, 1);}
  return Rcpp::List::create(
    _["x"] = x, _["y"] = y, _["rad"] = r,
    _["centroid"] = centroid, _["clRad"] = clrad
  );
}

// function to abstract the overlap checking and refitting process if an initial
// circle fitting attempt overlaps during iteration
int fit_circle(
  Node& curr_circ, Node& nxt_circ, std::vector<Node>& circles, int& j,
  int interupt_check_modulo, bool progbar, int num_circles
) {
  std::pair<Node*, Node*> check = overlap_check(
    curr_circ, nxt_circ, circles[j]
  );
  
  if (is_clear(check)) {
    insert_circle(curr_circ, nxt_circ, circles[j]);
    j++;
    periodic_interupt_check(j, interupt_check_modulo);
    if (progbar && (j <= num_circles)) {
      progress_bar(j, num_circles);
    }
  } else {
    while (!is_clear(check)) {
      Node& Cm = *(check.first);
      Node& Cn = *(check.second);
      
      fwd_remove(Cm, Cn);
      fit_tang_circle(Cm, Cn, circles[j]);
      check = overlap_check(Cm, Cn, circles[j]);
      
      if (is_clear(check)) {
        insert_circle(Cm, Cn, circles[j]);
        j++;
        periodic_interupt_check(j, interupt_check_modulo);
      }
    }
  }
  return j;
}

// this can be removed soon as its been accounted for in the subsequent function
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

// return the R list from a packed vector of circles for at least 3 circles
Rcpp::List process_into_clusterlist(
    std::vector<Node>& circles,
    Rcpp::NumericVector& centroid,
    double rad_scale_factor,
    int num_circles,
    bool progbar
) {
  Rcpp::NumericVector x (num_circles), y (num_circles), rad (num_circles);
  bool do_shift_x = (centroid[0] != 0), do_shift_y = (centroid[1] != 0);
  bool do_scale_rad = (rad_scale_factor != 1);
  int max_x_index = 0;
  
  for (int i = 0; i < num_circles; i++) {
    if (do_shift_x) {
      x[i] = circles[i].x + centroid[0];
    } else {
      x[i] = circles[i].x;
    }
    if (x[i] > x[max_x_index]) {
      max_x_index = i;
    }
    if (do_shift_y) {
      y[i] = circles[i].y + centroid[1];
    } else {
      y[i] = circles[i].y;
    }
    if (do_scale_rad) {
      rad[i] = circles[i].rad * rad_scale_factor;
    } else {
      rad[i] = circles[i].rad;
    }
  }
  double clRad = x[max_x_index] + (rad[max_x_index]/rad_scale_factor);
  if (progbar) {progress_bar(1, 1);}
  return Rcpp::List::create(
    _["x"] = x, _["y"] = y, _["rad"] = rad, _["centroid"] = centroid,
    _["clRad"] = clRad
  );
}

// Not exported because it causes an infinite loop atm :/
Rcpp::List circle_layout(
  std::vector<double> input_rad_vec,
  Rcpp::NumericVector centroid,
  double rad_scale_factor = 1,
  bool ORDER = true,
  bool try_place = false,
  bool progbar = true
) {
  if (progbar) {progress_bar(0, 1);}
  int num_circles = input_rad_vec.size();
  
  // if ORDER, decreasing sort of input_rad_vec
  if (ORDER && (num_circles > 1)) {
    std::sort(all(input_rad_vec), std::greater<int>());
  }
  
  // initialize list of circle nodes with radii
  std::vector<Node> circles (num_circles);
  for (int i = 0; i < num_circles; i++) {
    circles[i] = Node(0, 0, input_rad_vec[i]);
  }
  
  if (is_degenerate_case(num_circles)) {
    return handle_degenerate_cases(
      num_circles, circles, centroid, rad_scale_factor, progbar
    );
  }
  
  // Place the first three circles to be mutually tangent,
  // with centroid at the origin, and link the nodes
  place_starting_three(circles[0], circles[1], circles[2]);
  init_boundary({&circles[0], &circles[1], &circles[2]});
  
  // loop through the remaining circles, fitting them
  int j = 3;
  int interupt_factor = int(std::round(num_circles * 0.1));
  
  while (j < num_circles) {
    Rcpp::Rcout << "\nloop " << j << "\n";
    
    periodic_interupt_check(j, interupt_factor);
    Node& curr_circ = circles[j];
    Rcpp::Rcout << "curr_circ done, preplacement\n";
    if (try_place) {
      curr_circ = closest_place(circles[j-1], circles[j]);
    } else {
      curr_circ = closest(circles[j-1]);
    }
    Rcpp::Rcout << "placed\n";
    
    // fit circle, check for overlap, and refit till condition satisfied
    Node& nxt_circ = *(curr_circ.nxt);
    fit_tang_circle(curr_circ, nxt_circ, circles[j]);
    
    // seems like the problem is actually here in overlap_check??? But sometimes not??
    Rcpp::Rcout << "pre-overlap_check done\n";

    // check for overlaps and update, refit, and recheck until clear.
    j = fit_circle(
      curr_circ, nxt_circ, circles, j, interupt_factor, progbar, num_circles
    );
  }
  return process_into_clusterlist(
    circles, centroid, rad_scale_factor, num_circles, progbar
  );
}

// I wonder whats the best practisce for these "intermediate functions". Maybe
// pass an input class around? But thats a bit clunky it feels like.