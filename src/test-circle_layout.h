#include <testthat.h>

// some of these tests should probably also test what happens to the pointers :/

// testing of two nodes' numeric parameters
bool xyr_are_equal(Node& c1, Node c2, double thr = 5e-5) {
  return approx_equal(c1.x, c2.x, thr) && 
    approx_equal(c1.y, c2.y, thr) && 
    approx_equal(c1.rad, c2.rad, thr);
}

// testdata, assuming init_boundary works
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

bool coords_are_equal(const Node& c1, const Node& c2) {
  return c1.x == c2.x && c1.y == c2.y && c1.rad == c2.rad;
}

context("Cpp circle_layout functions") {
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
    expect_true(closest(nodes.c3) == &nodes.c1);
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
  
  // unfinished
  test_that("overlap_check works") {
    expect_true(1==1);
  }
  
  test_that("is_degenerate_case works") {
    expect_true(is_degenerate_case(1));
    expect_true(is_degenerate_case(2));
    expect_false(is_degenerate_case(3));
  }
  
  // unfinished
  test_that("handle_degenerate_cases works") {
    expect_true(1==1);
  }
}
