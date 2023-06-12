// function for testing of two nodes' numeric parameters x, y, rad.
bool xyr_are_equal(Node& c1, Node c2, double thr = 5e-5) {
  return approx_equal(c1.x, c2.x, thr) && 
    approx_equal(c1.y, c2.y, thr) && 
    approx_equal(c1.rad, c2.rad, thr);
}

// some linked test nodes for testing, assuming init_boundary works
class Testdata {
  public:
    Node c1; Node c2; Node c3; Node c4; Node c5;
    Testdata(
      const Node& a = Node(1, 2, 3),
      const Node& b = Node(6, 7, 8),
      const Node& c = Node(50, 90, 1),
      const Node& d = Node(-1, -3, 3),
      const Node& e = Node(-10, -69, 4)
    ) : c1(a), c2(b), c3(c), c4(d), c5(e) {
      init_boundary({&c1, &c2, &c3, &c4, &c5});
    }
};

context("C++ circle_layout") {
  test_that("init_boundary works") {
    Node c1 = Node(1, 2, 3);
    Node c2 = Node(6, 7, 8);
    Node c3 = Node(50, 90, 1);
    std::vector<Node*> nodes = {&c1, &c2, &c3};
    
    init_boundary(nodes);
    expect_true(c1.nxt == &c2); expect_true(c1.prv == &c3);
    expect_true(c2.nxt == &c3); expect_true(c2.prv == &c1);
    expect_true(c3.nxt == &c1); expect_true(c3.prv == &c2);
  }

  test_that("insert_circle_works") {
    Node c1 = Node(1, 2, 3);
    Node c3 = Node(6, 7, 8);
    c1.nxt = &c3; c1.prv = &c3;
    c3.nxt = &c1, c3.prv = &c1;
    
    Node c2 = Node(50, 90, 1);
    insert_circle(c1, c3, c2);
    
    expect_true(c1.nxt == &c2);
    expect_true(c3.prv == &c2);
    expect_true(c2.nxt == &c3);
    expect_true(c2.prv == &c1);
  }
  
  test_that("fwd_dist works") {
    Node c1 = Node(1, 2, 3);
    Node c2 = Node(6, 7, 8);
    Node c3 = Node(50, 90, 1);
    init_boundary({&c1, &c2, &c3});
    
    expect_true(fwd_dist(c3, c1) == 1);
    expect_true(fwd_dist(c1, c2) == 1);
    expect_true(fwd_dist(c2, c1) == 2);
  }
 
  test_that("fwd_remove works") {
    Testdata nodes = Testdata();
   
    fwd_remove(nodes.c2, nodes.c5);
    expect_true(nodes.c2.nxt == &nodes.c5);
    expect_true(nodes.c5.prv == &nodes.c2);
    expect_true(nodes.c5.nxt == &nodes.c1);
  }
  
  test_that("centre_dist works") {
    Node c1 = Node(3, 4, 1);
    expect_true(centre_dist(c1) == 5);
    
    Node c2 = Node(-12, -35, 1);
    expect_true(centre_dist(c2) == 37);
  }
  
  test_that("fit_tang_circle works") {
    Node c1 = Node(1, 2, 3);
    Node c2 = Node(6, 7, 8);
    Node c3 = Node(50, 90, 1);
    init_boundary({&c1, &c2, &c3});
    
    expect_true(&fit_tang_circle(c1, c2, c3) == &c3);
    expect_true(approx_equal(c3.x, -2.47718));
    expect_true(approx_equal(c3.y, 3.97718));
  }
  
  test_that("tang_circle_dist works") {
    Node c1 = Node(1, 2, 3);
    Node c2 = Node(6, 7, 8);
    Node c3 = Node(50, 90, 1);
    init_boundary({&c1, &c2, &c3});
    
    expect_true(approx_equal(tang_circle_dist(c1, c2, c3), 4.685548));
  }
  
  test_that("place_starting_three works") {
    Node c1 = Node(1, 2, 3);
    Node c2 = Node(6, 7, 8);
    Node c3 = Node(50, 90, 1);
    init_boundary({&c1, &c2, &c3});
    
    place_starting_three(c1, c2, c3);
    expect_true(xyr_are_equal(c1, Node(-4.984898, -1.466558, 3)));
    expect_true(xyr_are_equal(c2, Node(6.015102, 3.533442, 8)));
    expect_true(xyr_are_equal(c3, Node(-1.030204, -2.066885, 1)));
  }
  
  test_that("closest works") {
    Testdata nodes = Testdata();
    expect_true(closest(nodes.c3) == nodes.c1);
    expect_true(closest(nodes.c1) == closest(nodes.c5));
    expect_true(closest(nodes.c2) == closest(nodes.c4));
  }
  
  // this is the suspicious function  - unfinished
  test_that("closest_place works") {
    expect_true(1==1);
  }
  
  test_that("do_intersect works") {
    Testdata nodes = Testdata();
    expect_true(do_intersect(nodes.c1, nodes.c2));
    expect_true(do_intersect(nodes.c1, nodes.c4));
    expect_false(do_intersect(nodes.c2, nodes.c5));
    expect_false(do_intersect(nodes.c3, nodes.c4));
  }
  
  test_that("geod_dist works") {
    Testdata nodes = Testdata();
    expect_true(geod_dist(nodes.c1, nodes.c2, nodes.c3) == 1);
    expect_true(geod_dist(nodes.c3, nodes.c2, nodes.c1) == 2);
    expect_true(geod_dist(nodes.c2, nodes.c5, nodes.c3) == 3);
    expect_true(geod_dist(nodes.c4, nodes.c1, nodes.c5) == 4);
  }
  
  test_that("make_sentinel_pair works") {
    std::pair<Node*, Node*> trial = make_sentinel_pair();
    expect_true(trial.first == nullptr);
    expect_true(trial.second == nullptr);
  }
  
  test_that("is_clear works") {
    std::pair<Node*, Node*> sentinel_pair = make_sentinel_pair();
    expect_true(is_clear(sentinel_pair));
    
    Node node1 = Node(1, 2, 3), node2 = Node(3, 2, 1);
    std::pair<Node*, Node*> false_pair = std::make_pair(&node1, &node2);
    expect_false(is_clear(false_pair));
  }
  
  test_that("overlap_check works for 5 connected nodes") {
    Testdata nodes = Testdata();
    std::pair<Node*, Node*> trial, expected;
    
    trial = overlap_check(nodes.c1, nodes.c2, nodes.c5);
    expected = std::make_pair(&nodes.c5, &nodes.c2);
    expect_true(trial == expected);
    
    trial = overlap_check(nodes.c3, nodes.c5, nodes.c4);
    expected = std::make_pair(&nodes.c3, &nodes.c1);
    expect_true(trial == expected);
    
    trial = overlap_check(nodes.c4, nodes.c2, nodes.c5);
    expect_true(is_clear(trial));
  }
  
  test_that("overlap_check works for 3 connected nodes") {
    Node n1 = Node(1,2,3), n2 = Node(6,7,8), n3 = Node(50,90,1);
    init_boundary({&n1, &n2, &n3});
    std::pair<Node*, Node*> trial, expected;
    trial = overlap_check(n1, n2, n3);
    
    expected = std::make_pair(&n1, &n3);
    expect_true(trial == expected);
  }
  
  test_that("is_degenerate_case works") {
    expect_true(is_degenerate_case(1));
    expect_true(is_degenerate_case(2));
    expect_false(is_degenerate_case(3));
  }
  
  test_that("handle_degenerate_cases for n = 1 works") {
    // generating test trial 1 values
    std::vector<Node> circles = {Node(0, 0, 69)};
    Rcpp::NumericVector centroid = {0, 0};
    Rcpp::List trial = handle_degenerate_cases(
      1, circles, centroid, 1
    );
    Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
    Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
    double trial_clRad = trial[4];
    
    // testing for equivalence to expected values for trial 1
    Rcpp::NumericVector x = {0}, y = {0}, rad = {69};

    expect_true(elements_are_equal(trial_x, x));
    expect_true(elements_are_equal(trial_y, y));
    expect_true(elements_are_equal(trial_rad, rad));
    expect_true(elements_are_equal(trial_centroid, centroid));
    expect_true(approx_equal(trial_clRad, 69.0));
    
    // generating test trial 2 values with new centroid and rad_scale_factor
    circles = {Node(0, 0, 69)};
    centroid = {4, 5};
    trial = handle_degenerate_cases(
      1, circles, centroid, 0.8
    );
    trial_x = trial[0], trial_y = trial[1], trial_rad = trial[2];
    trial_centroid = trial[3], trial_clRad = trial[4];
    
    // testing for equivalence to expected values for trial 1
    x = {4}, y = {5}, rad = {69 * 0.8};
    
    expect_true(elements_are_equal(trial_x, x));
    expect_true(elements_are_equal(trial_y, y));
    expect_true(elements_are_equal(trial_rad, rad));
    expect_true(elements_are_equal(trial_centroid, centroid));
    expect_true(approx_equal(trial_clRad, 69 * 0.8));
  }
  
  test_that("handle_degenerate_cases for n = 2 works") {
    // generating test trial 1 values
    std::vector<Node> circles = {Node(0, 0, 69), Node(0, 0, 420)};
    Rcpp::NumericVector centroid = {0, 0};
    Rcpp::List trial = handle_degenerate_cases(
      2, circles, centroid, 1
    );
    Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
    Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
    double trial_clRad = trial[4];
    
    // testing for equivalence to expected values for trial 1
    Rcpp::NumericVector x = {-69, 420}, y = {0, 0}, rad = {69, 420};
    
    expect_true(elements_are_equal(trial_x, x));
    expect_true(elements_are_equal(trial_y, y));
    expect_true(elements_are_equal(trial_rad, rad));
    expect_true(elements_are_equal(trial_centroid, centroid));
    expect_true(approx_equal(trial_clRad, 244.5));
    
    // generating test trial 2 values with new centroid and rad_scale_factor
    circles = {Node(0, 0, 4), Node(0, 0, 3)};
    centroid = {1, 2};
    trial = handle_degenerate_cases(
      2, circles, centroid, 0.8
    );
    trial_x = trial[0], trial_y = trial[1], trial_rad = trial[2];
    trial_centroid = trial[3], trial_clRad = trial[4];
    
    // testing for equivalence to expected values for trial 1
    x = {-3, 4}, y = {2, 2}, rad = {3.2, 2.4};
    
    expect_true(elements_are_equal(trial_x, x));
    expect_true(elements_are_equal(trial_y, y));
    expect_true(elements_are_equal(trial_rad, rad));
    expect_true(elements_are_equal(trial_centroid, centroid));
    expect_true(approx_equal(trial_clRad, 2.8));
  }

  test_that("process_into_clusterlist works") {
    // create trial data, 5 circles, centroid (1,1), rad_scale_factor 0.8
    std::vector<Node> circles = {
      Node(1,2,3), Node(6,7,8), Node(50,90,1), Node(-1,-3,3), Node(-10,-69,4)
    };
    Rcpp::NumericVector centroid = {1, 1};
    Rcpp::List trial = process_into_clusterlist(
      circles, centroid, 0.8, 5, false
    );
    Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
    Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
    double trial_clRad = trial[4];
    
    // test for equivalence to expected data
    Rcpp::NumericVector x = {2, 7, 51, 0, -9};
    Rcpp::NumericVector y = {3, 8, 91, -2, -68};
    Rcpp::NumericVector rad = {2.4, 6.4, 0.8, 2.4, 3.2};
    
    expect_true(elements_are_equal(trial_x, x));
    expect_true(elements_are_equal(trial_y, y));
    expect_true(elements_are_equal(trial_rad, rad));
    expect_true(elements_are_equal(trial_centroid, centroid));
    expect_true(approx_equal(trial_clRad, 51));
  }
  
  test_that("circle_layout for 3 circles works") {
    Rcpp::NumericVector centroid = {4, 3};
    Rcpp::List trial = circle_layout({6,9,2}, centroid, 1, true, false, false);
    Rcpp::NumericVector trial_x = trial[0], trial_y = trial[1];
    Rcpp::NumericVector trial_rad = trial[2], trial_centroid = trial[3];
    double trial_clRad = trial[4];
    
    Rcpp::NumericVector x = {-4.133333, 10.866667, 5.266667};
    Rcpp::NumericVector y = {4.9043809, 4.9043809, -0.8087618};
    Rcpp::NumericVector rad = {9, 6, 2};
    
    expect_true(elements_are_equal(trial_x, x));
    expect_true(elements_are_equal(trial_y, y));
    expect_true(elements_are_equal(trial_rad, rad));
    expect_true(elements_are_equal(trial_centroid, centroid));
    expect_true(approx_equal(trial_clRad, 12.86667));
  }
  
  // circle layout needs to be tested more to avoid the infinite loop. Could also be done in R
}

// some of these tests should probably also test what happens to the pointers :/