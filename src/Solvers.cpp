#include "truss/Solvers.h"

#include <cmath>
#include <sstream>

namespace truss {

namespace {

// Below this size the OpenMP thread start-up costs more than the loop itself.
// Only the '#pragma omp' lines use it, hence [[maybe_unused]] for serial builds.
[[maybe_unused]] constexpr int kParallelThreshold = 256;

// A pivot smaller than this fraction of the original diagonal entry is treated as zero.
constexpr double kPivotTolerance = 1e-10;

int size_of(const std::vector<double>& v) { return static_cast<int>(v.size()); }

}  // namespace

std::string to_string(SolverType type) {
    return type == SolverType::Gauss ? "gauss" : "cg";
}

std::vector<double> multiply(const Matrix& A, const std::vector<double>& x) {
    const int n = size_of(x);
    std::vector<double> y(x.size(), 0.0);
#pragma omp parallel for schedule(static) if (n > kParallelThreshold)
    for (int i = 0; i < n; ++i) {
        const std::vector<double>& row = A[i];
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += row[j] * x[j];
        }
        y[i] = sum;
    }
    return y;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    const int n = size_of(a);
    double sum = 0.0;
#pragma omp parallel for reduction(+ : sum) schedule(static) if (n > kParallelThreshold)
    for (int i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<double> gauss_solve(Matrix A, std::vector<double> b) {
    const int n = size_of(b);

    std::vector<double> scale(b.size());
    for (int i = 0; i < n; ++i) {
        scale[i] = std::abs(A[i][i]) > 0.0 ? std::abs(A[i][i]) : 1.0;
    }

    // Forward elimination to upper triangular form.
    for (int k = 0; k < n; ++k) {
        const double pivot = A[k][k];
        if (!(std::abs(pivot) > kPivotTolerance * scale[k])) {
            throw SingularMatrixError(k);
        }
        const std::vector<double>& row_k = A[k];

        // Rows below k are independent of each other, so they are updated in parallel.
        // Truss matrices are banded: most rows have a zero in column k and are skipped,
        // which is why the schedule is dynamic.
#pragma omp parallel for schedule(dynamic, 16) if (n - k > kParallelThreshold)
        for (int i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / pivot;
            if (factor == 0.0) {
                continue;
            }
            std::vector<double>& row_i = A[i];
            for (int j = k + 1; j < n; ++j) {
                row_i[j] -= factor * row_k[j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution.
    std::vector<double> x(b.size(), 0.0);
    for (int k = n - 1; k >= 0; --k) {
        double sum = b[k];
        for (int j = k + 1; j < n; ++j) {
            sum -= A[k][j] * x[j];
        }
        x[k] = sum / A[k][k];
    }
    return x;
}

std::vector<double> cg_solve(const Matrix& A, const std::vector<double>& b, double tol,
                             int max_iter, int* iterations) {
    const int n = size_of(b);
    std::vector<double> x(b.size(), 0.0);
    std::vector<double> r = b;
    std::vector<double> p = r;

    double rr = dot(r, r);
    const double b_norm = std::sqrt(rr);
    if (iterations != nullptr) {
        *iterations = 0;
    }
    if (b_norm == 0.0) {
        return x;  // no load: zero displacement
    }

    for (int iter = 1; iter <= max_iter; ++iter) {
        const std::vector<double> Ap = multiply(A, p);
        const double pAp = dot(p, Ap);
        if (!(pAp > 0.0)) {
            throw std::runtime_error(
                "CG breakdown: the stiffness matrix is not positive definite "
                "(the truss is probably a mechanism)");
        }
        const double alpha = rr / pAp;
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        const double rr_new = dot(r, r);
        if (std::sqrt(rr_new) <= tol * b_norm) {
            if (iterations != nullptr) {
                *iterations = iter;
            }
            return x;
        }
        const double beta = rr_new / rr;
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }
        rr = rr_new;
    }

    std::ostringstream msg;
    msg << "CG did not converge in " << max_iter << " iterations (relative residual "
        << std::sqrt(rr) / b_norm << "); the truss may be a mechanism";
    throw std::runtime_error(msg.str());
}

}  // namespace truss
