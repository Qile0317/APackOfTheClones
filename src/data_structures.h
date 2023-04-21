#include <vector>
#include <string>
#include <iostream>

using namespace std;

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

    // printer and traverser for debugging
    friend ostream& operator<<(std::ostream& os, const Node& node) {
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

// a cluster of circles
class Cluster {
public:
  vector<double> x;             
  vector<double> y;             
  vector<double> rad;           
  vector<double> centroid;      
  double cluster_rad;                 
  string color;                
  
  Cluster(vector<double> x_val, vector<double> y_val, vector<double> rad_val,
          vector<double> centroid_val, double clRad_val, string color_val) {
    x = x_val;
    y = y_val;
    rad = rad_val;
    centroid = centroid_val;
    cluster_rad = clRad_val;
    color = color_val;
  }
};