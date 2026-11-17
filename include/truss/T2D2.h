#pragma once

#include <array>
#include <vector>

#include "truss/Node.h"

namespace truss {

using Block2 = std::array<std::array<double, 2>, 2>;
using Matrix4 = std::array<std::array<double, 4>, 4>;

// Two-node linear truss element in 2D (named after the Abaqus element T2D2).
//
// With c = cos(theta), s = sin(theta) of the axis from node 1 to node 2:
//
//          EA | c*c  c*s |                    |  k  -k |
//     k =  -- |          |       K_element =  |        |
//           L | c*s  s*s |                    | -k   k |
//
// Length and direction are computed once in the constructor, so the element is
// read-only afterwards and safe to use from several threads at the same time.
class T2D2 {
public:
    T2D2(int id, const Node& node1, const Node& node2, double E, double A);

    int id() const { return id_; }
    int node1() const { return node1_; }
    int node2() const { return node2_; }
    double E() const { return E_; }
    double A() const { return A_; }
    double length() const { return length_; }
    double cos() const { return c_; }
    double sin() const { return s_; }

    // The 2x2 block k above.
    Block2 stiffness_block() const;

    // The full 4x4 element stiffness for the DOFs (u1x, u1y, u2x, u2y).
    Matrix4 stiffness() const;

    // Axial force N = EA/L * elongation, from the global displacement vector.
    // Positive in tension, negative in compression.
    double axial_force(const std::vector<double>& u) const;

    double stress(const std::vector<double>& u) const { return axial_force(u) / A_; }

private:
    int id_;
    int node1_;
    int node2_;
    double E_;
    double A_;
    double length_;
    double c_;
    double s_;
};

}  // namespace truss
