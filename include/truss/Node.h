#pragma once

#include <vector>

namespace truss {

// A pin joint of the truss. Node i owns the global DOFs 2*i (ux) and 2*i + 1 (uy).
class Node {
public:
    Node(int id, double x, double y) : id_{id}, x_{x}, y_{y} {}

    int id() const { return id_; }
    double x() const { return x_; }
    double y() const { return y_; }

    // Elements connected to this node. Assembly runs node by node over this list.
    const std::vector<int>& elements() const { return elements_; }
    void add_element(int element_id) { elements_.push_back(element_id); }

private:
    int id_;
    double x_;
    double y_;
    std::vector<int> elements_;
};

}  // namespace truss
