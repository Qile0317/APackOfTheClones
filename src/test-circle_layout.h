#include <testthat.h>

class Testdata {
public:
  Node c1;
  Node c2;
  Node c3;
  
  Testdata(
    const Node& a = Node(1, 2, 3),
    const Node& b = Node(6, 7, 8),
    const Node& c = Node(50, 90, 1)
  ) : c1(a), c2(b), c3(c) {
    c1.nxt = &c2; c1.prv = &c3;
    c2.nxt = &c3; c2.prv = &c1;
    c3.nxt = &c1; c3.prv = &c2;
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
    expect_true(c1.nxt == &c2);
    expect_true(c1.prv == &c3);
    expect_true(c2.nxt == &c3);
    expect_true(c2.prv == &c1);
    expect_true(c3.nxt == &c1);
    expect_true(c3.prv == &c2);
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
    Testdata foo = Testdata();
    expect_true(fwd_dist(foo.c3, foo.c1) == 1);
    expect_true(fwd_dist(foo.c1, foo.c2) == 1);
    expect_true(fwd_dist(foo.c2, foo.c1) == 2);
  }
 
  test_that("fwd_remove works") {
    Node c1 = Node(1, 2, 3);
    Node c2 = Node(6, 7, 8);
    Node c3 = Node(50, 90, 1);
    Node c4 = Node(-1, -2, -3);
    Node c5 = Node(-10, -69, 4);
    std::vector<Node*> nodes = {&c1, &c2, &c3, &c4, &c5};
    init_boundary(nodes);
   
    fwd_remove(c2, c5);
    expect_true(c2.nxt == &c5);
    expect_true(c5.prv == &c2);
    expect_true(c5.nxt == &c1);
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
    expect_true(approx_equal(c3.x, -2.47718, 5e-5));
    expect_true(approx_equal(c3.y, 3.97718, 5e-5));
  }
}
