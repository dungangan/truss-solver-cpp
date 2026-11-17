#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "truss/Node.h"
#include "truss/Solvers.h"
#include "truss/T2D2.h"

namespace truss {

// A 2D truss: nodes, T2D2 elements, nodal loads and pinned/roller supports.
//
// Typical use:
//     Model model = Model::from_file("truss.txt");
//     model.assemble();
//     model.solve(SolverType::Gauss);
//     model.print_results(std::cout);
class Model {
public:
    // Read a model (format in README.md). Errors report the file name and line.
    static Model from_file(const std::string& path);
    static Model from_stream(std::istream& in, const std::string& source_name = "<input>");

    // Assemble the global stiffness matrix K in parallel with OpenMP.
    void assemble();

    // Apply the supports, solve K u = f, then compute member forces and reactions.
    // Calls assemble() first if needed.
    SolveInfo solve(SolverType type);

    int num_nodes() const { return static_cast<int>(nodes_.size()); }
    int num_elements() const { return static_cast<int>(elements_.size()); }
    int num_dofs() const { return 2 * num_nodes(); }
    int num_free_dofs() const;

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<T2D2>& elements() const { return elements_; }
    bool is_fixed(int dof) const { return fixed_[dof] != 0; }

    const std::vector<double>& loads() const { return f_; }
    const Matrix& stiffness() const { return K_; }  // assembled K, before supports

    // Results, available after solve().
    const std::vector<double>& displacements() const { return u_; }
    const std::vector<double>& axial_forces() const { return axial_forces_; }
    const std::vector<double>& reactions() const { return reactions_; }  // zero at free DOFs

    // Sum of all external forces (loads + reactions) in x and y. Close to zero for a
    // correct solution.
    double force_balance_x() const;
    double force_balance_y() const;

    void print_input(std::ostream& os) const;
    void print_results(std::ostream& os) const;
    void print_matrices(std::ostream& os) const;

private:
    Model() = default;
    void compute_results();

    std::vector<Node> nodes_;
    std::vector<T2D2> elements_;
    std::vector<double> f_;     // nodal loads, size 2n
    std::vector<char> fixed_;   // 1 if the DOF is supported, size 2n
    Matrix K_;                  // global stiffness, size 2n x 2n

    bool solved_ = false;
    std::vector<double> u_;
    std::vector<double> axial_forces_;
    std::vector<double> reactions_;
};

}  // namespace truss
