#include "circle_packing.h"

// for testing
int main() {
  Node c1 = Node(1,2,3);
  Node c3 = Node(6,7,8);
  Node c4 = Node(50,90,1);
  init_boundary(c1,c3,c4);
  cout << fwd_dist(c4,c3) << endl; // expect 2

  Node c2 = Node(100,100,3);
  insert_circle(c1,c3,c2);
  cout << fwd_dist(c2,c1) << endl; // expect 3
  cout << c1 << endl;
  cout << "-----------" << endl;

  Node c5 = Node(20,200,3.5);
  insert_circle(c4,c1,c5);
  c1.traverse();
  cout << "-----------" << endl;

  fwd_remove(c3,c5);
  c1.traverse();
  cout << "-----------" << endl;

  cout << centre_dist(c1) << endl;
  cout << "-----------" << endl;

  c1 = Node(0,0,1); c2 = Node(2,0,1); c3 = Node(0,2,1);
  init_boundary(c1,c3,c4); insert_circle(c1,c3,c2); insert_circle(c4,c1,c5);
  c1.traverse();
  cout << "after" << endl;
  place_starting_three(c1,c2,c3);
  c1.traverse();
  cout << "-----------" << endl;
  return 0;
}