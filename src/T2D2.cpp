#include "truss/T2D2.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace truss {

T2D2::T2D2(int id, const Node& node1, const Node& node2, double E, double A)
    : id_{id}, node1_{node1.id()}, node2_{node2.id()}, E_{E}, A_{A} {
    if (!(E > 0.0) || !(A > 0.0)) {
        throw std::invalid_argument("element " + std::to_string(id) +
                                    ": E and A must be positive");
    }
    const double dx = node2.x() - node1.x();
    const double dy = node2.y() - node1.y();
    length_ = std::hypot(dx, dy);
    if (!(length_ > 0.0)) {
        throw std::invalid_argument("element " + std::to_string(id) +
                                    " has zero length (its nodes are at the same position)");
    }
    c_ = dx / length_;
    s_ = dy / length_;
}

Block2 T2D2::stiffness_block() const {
    const double k = E_ * A_ / length_;
    return {{{k * c_ * c_, k * c_ * s_},
             {k * c_ * s_, k * s_ * s_}}};
}

Matrix4 T2D2::stiffness() const {
    const Block2 k = stiffness_block();
    Matrix4 K{};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            K[i][j] = k[i][j];
            K[i + 2][j + 2] = k[i][j];
            K[i][j + 2] = -k[i][j];
            K[i + 2][j] = -k[i][j];
        }
    }
    return K;
}

double T2D2::axial_force(const std::vector<double>& u) const {
    const double dux = u[2 * node2_] - u[2 * node1_];
    const double duy = u[2 * node2_ + 1] - u[2 * node1_ + 1];
    const double elongation = c_ * dux + s_ * duy;
    return E_ * A_ / length_ * elongation;
}

}  // namespace truss
