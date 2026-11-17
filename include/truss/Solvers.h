#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace truss {

// Dense row-major matrix. Fine for trusses up to a few thousand nodes.
using Matrix = std::vector<std::vector<double>>;

enum class SolverType { Gauss, CG };

std::string to_string(SolverType type);

struct SolveInfo {
    SolverType type = SolverType::Gauss;
    int iterations = 0;    // CG iterations (0 for Gauss)
    double residual = 0.0; // relative residual |K u - f| / |f| of the solved system
};

// Thrown when elimination meets a zero pivot. `row` is the DOF where it happened.
class SingularMatrixError : public std::runtime_error {
public:
    explicit SingularMatrixError(int row)
        : std::runtime_error("matrix is singular at row " + std::to_string(row)), row_{row} {}
    int row() const { return row_; }

private:
    int row_;
};

// Gaussian elimination without pivoting. The stiffness matrix is symmetric positive
// definite once the supports are applied, so pivoting is not needed; a zero pivot
// means the truss is a mechanism and raises SingularMatrixError.
std::vector<double> gauss_solve(Matrix A, std::vector<double> b);

// Conjugate gradient method for symmetric positive definite A, starting from x = 0.
// Stops when |r| <= tol * |b|. Throws std::runtime_error if it breaks down or does
// not converge within max_iter iterations.
std::vector<double> cg_solve(const Matrix& A, const std::vector<double>& b, double tol,
                             int max_iter, int* iterations = nullptr);

// Helpers shared by the solvers and the model (parallel with OpenMP for large n).
std::vector<double> multiply(const Matrix& A, const std::vector<double>& x);
double dot(const std::vector<double>& a, const std::vector<double>& b);

}  // namespace truss
