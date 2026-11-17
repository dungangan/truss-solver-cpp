// Tests for the truss solver. Plain C++ with no test framework: each test is a
// function, CHECK macros count failures, and the exit code is non-zero if any failed.

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "truss/Model.h"

using truss::Model;
using truss::SolverType;

namespace {

int g_checks = 0;
int g_failures = 0;

void report_failure(const char* file, int line, const std::string& message) {
    ++g_failures;
    std::cerr << file << ":" << line << ": FAILED: " << message << '\n';
}

void check_close(double actual, double expected, double tolerance, const char* expr,
                 const char* file, int line) {
    ++g_checks;
    const double scale = std::max(1.0, std::abs(expected));
    if (!(std::abs(actual - expected) <= tolerance * scale)) {
        std::ostringstream msg;
        msg.precision(15);
        msg << expr << " = " << actual << ", expected " << expected;
        report_failure(file, line, msg.str());
    }
}

#define CHECK(cond)                                                  \
    do {                                                             \
        ++g_checks;                                                  \
        if (!(cond)) report_failure(__FILE__, __LINE__, #cond);      \
    } while (0)

// Relative tolerance for values larger than 1, absolute below.
#define CHECK_CLOSE(actual, expected, tol) \
    check_close((actual), (expected), (tol), #actual, __FILE__, __LINE__)

// Checks that `expr` throws an exception whose message contains `text`.
#define CHECK_THROWS_WITH(expr, text)                                                 \
    do {                                                                              \
        ++g_checks;                                                                   \
        try {                                                                         \
            (void)(expr);                                                             \
            report_failure(__FILE__, __LINE__, std::string(#expr) + " did not throw"); \
        } catch (const std::exception& e) {                                           \
            if (std::string(e.what()).find(text) == std::string::npos) {              \
                report_failure(__FILE__, __LINE__,                                    \
                               std::string("message '") + e.what() +                  \
                                   "' does not contain '" + (text) + "'");            \
            }                                                                         \
        }                                                                             \
    } while (0)

Model parse(const std::string& text) {
    std::istringstream in(text);
    return Model::from_stream(in, "test");
}

const char* const kFourNodeTruss = R"(
4 5 1 2
0 0 0
1 1000 0
2 0 1500
3 1000 1500
0 0 2 200000 225
1 2 3 200000 225
2 0 3 200000 225
3 2 1 200000 225
4 1 3 200000 225
2 50000 0
0 1 1
1 1 1
)";

// ---------------------------------------------------------------------------

void test_single_bar() {
    // Horizontal bar, pinned at node 0, roller (fixed in y) at node 1, pulled by F.
    // Exact: u = F L / (E A), N = F, reaction = -F.
    const double L = 2000.0, E = 210000.0, A = 100.0, F = 10000.0;
    std::ostringstream text;
    text << "2 1 1 2\n0 0 0\n1 " << L << " 0\n0 0 1 " << E << " " << A << "\n1 " << F
         << " 0\n0 1 1\n1 0 1\n";
    for (const SolverType solver : {SolverType::Gauss, SolverType::CG}) {
        Model model = parse(text.str());
        model.solve(solver);
        CHECK_CLOSE(model.displacements()[2], F * L / (E * A), 1e-12);
        CHECK_CLOSE(model.displacements()[3], 0.0, 1e-12);
        CHECK_CLOSE(model.axial_forces()[0], F, 1e-12);
        CHECK_CLOSE(model.reactions()[0], -F, 1e-12);
    }
}

void test_two_bar_truss() {
    // Symmetric two-bar truss hanging from two pins, vertical load P at the apex.
    // Bars at angle theta with sin(theta) = 0.8, L = 2500:
    //   N = P / (2 sin(theta)),  deflection = P L / (2 E A sin^2(theta)).
    const double P = 40000.0, E = 200000.0, A = 500.0, L = 2500.0, s = 0.8;
    const std::string text = R"(
3 2 1 2
0 -1500 0
1  1500 0
2     0 -2000
0 0 2 200000 500
1 1 2 200000 500
2 0 -40000
0 1 1
1 1 1
)";
    for (const SolverType solver : {SolverType::Gauss, SolverType::CG}) {
        Model model = parse(text);
        model.solve(solver);
        CHECK_CLOSE(model.displacements()[4], 0.0, 1e-12);
        CHECK_CLOSE(model.displacements()[5], -P * L / (2 * E * A * s * s), 1e-12);
        CHECK_CLOSE(model.axial_forces()[0], P / (2 * s), 1e-12);
        CHECK_CLOSE(model.axial_forces()[1], P / (2 * s), 1e-12);
        CHECK_CLOSE(model.reactions()[0], -15000.0, 1e-12);
        CHECK_CLOSE(model.reactions()[1], 20000.0, 1e-12);
        CHECK_CLOSE(model.reactions()[2], 15000.0, 1e-12);
        CHECK_CLOSE(model.reactions()[3], 20000.0, 1e-12);
    }
}

void test_four_node_truss_reference() {
    // Statically indeterminate example (examples/four_node_truss.txt). Reference
    // values from an independent NumPy solution.
    const double u_ref[] = {5.393520947324e+00, 1.314207793613e+00, 4.866502188930e+00,
                            -1.185792206387e+00};
    const double N_ref[] = {3.942623380838e+04, -2.371584412775e+04, 4.275434602175e+04,
                            -4.738443586485e+04, -3.557376619162e+04};
    const double R_ref[] = {-2.371584412775e+04, -7.500000000000e+04, -2.628415587225e+04,
                            7.500000000000e+04};

    for (const SolverType solver : {SolverType::Gauss, SolverType::CG}) {
        Model model = parse(kFourNodeTruss);
        const truss::SolveInfo info = model.solve(solver);
        CHECK(info.residual < 1e-10);
        for (int d = 0; d < 4; ++d) CHECK_CLOSE(model.displacements()[d], 0.0, 1e-14);
        for (int d = 0; d < 4; ++d) CHECK_CLOSE(model.displacements()[4 + d], u_ref[d], 1e-10);
        for (int e = 0; e < 5; ++e) CHECK_CLOSE(model.axial_forces()[e], N_ref[e], 1e-10);
        for (int d = 0; d < 4; ++d) CHECK_CLOSE(model.reactions()[d], R_ref[d], 1e-10);
        CHECK_CLOSE(model.force_balance_x(), 0.0, 1e-8);
        CHECK_CLOSE(model.force_balance_y(), 0.0, 1e-8);
        if (solver == SolverType::CG) {
            CHECK(info.iterations > 0 && info.iterations <= model.num_free_dofs() + 2);
        }
    }
}

void test_element_stiffness() {
    // 30 degree bar: symmetric, no force under rigid translation, and EA/L for a
    // unit stretch along the axis.
    const double E = 70000.0, A = 50.0, L = 800.0;
    const double pi = std::acos(-1.0);
    const double c = std::cos(pi / 6.0), s = std::sin(pi / 6.0);
    const truss::Node n1(0, 10.0, 20.0);
    const truss::Node n2(1, 10.0 + L * c, 20.0 + L * s);
    const truss::T2D2 element(0, n1, n2, E, A);
    const truss::Matrix4 K = element.stiffness();

    CHECK_CLOSE(element.length(), L, 1e-12);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) CHECK_CLOSE(K[i][j], K[j][i], 1e-12);

    const double rigid[2][4] = {{1, 0, 1, 0}, {0, 1, 0, 1}};
    for (const auto& u : rigid) {
        for (int i = 0; i < 4; ++i) {
            double force = 0.0;
            for (int j = 0; j < 4; ++j) force += K[i][j] * u[j];
            CHECK_CLOSE(force, 0.0, 1e-9);
        }
    }

    // Stretch by 1: node 2 moves along the axis.
    const std::vector<double> u = {0.0, 0.0, c, s};
    CHECK_CLOSE(element.axial_force(u), E * A / L, 1e-12);
}

void test_assembly_symmetric_and_balanced() {
    Model model = parse(kFourNodeTruss);
    model.assemble();
    const truss::Matrix& K = model.stiffness();
    const int n = model.num_dofs();
    for (int i = 0; i < n; ++i) {
        double row_sum_x = 0.0, row_sum_y = 0.0;  // forces for rigid translation in x / y
        for (int j = 0; j < n; ++j) {
            CHECK_CLOSE(K[i][j], K[j][i], 1e-12);
            (j % 2 == 0 ? row_sum_x : row_sum_y) += K[i][j];
        }
        CHECK_CLOSE(row_sum_x, 0.0, 1e-9);
        CHECK_CLOSE(row_sum_y, 0.0, 1e-9);
    }
}

void test_large_grid_parallel_paths() {
    // X-braced grid, 15 x 10 nodes = 300 DOFs: large enough to run the OpenMP
    // branches of assembly, Gauss elimination and CG.
    const int nx = 15, ny = 10;
    const double h = 1000.0;
    auto id = [&](int i, int j) { return j * nx + i; };
    std::vector<std::string> elements;
    auto add = [&](int a, int b) {
        elements.push_back(std::to_string(elements.size()) + " " + std::to_string(a) + " " +
                           std::to_string(b) + " 200000 " + std::to_string(100 + a % 7));
    };
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (i + 1 < nx) add(id(i, j), id(i + 1, j));
            if (j + 1 < ny) add(id(i, j), id(i, j + 1));
            if (i + 1 < nx && j + 1 < ny) {
                add(id(i, j), id(i + 1, j + 1));
                add(id(i + 1, j), id(i, j + 1));
            }
        }
    }
    std::ostringstream text;
    text << nx * ny << " " << elements.size() << " " << nx << " " << nx << "\n";
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) text << id(i, j) << " " << i * h << " " << j * h << "\n";
    for (const std::string& e : elements) text << e << "\n";
    for (int i = 0; i < nx; ++i) text << id(i, ny - 1) << " " << 1000.0 * i << " -5000\n";
    for (int i = 0; i < nx; ++i) text << id(i, 0) << " 1 1\n";

    // The parallel row-wise assembly must equal the textbook element-by-element scatter.
    Model model = parse(text.str());
    model.assemble();
    const int n = model.num_dofs();
    truss::Matrix K_ref(n, std::vector<double>(n, 0.0));
    for (const truss::T2D2& e : model.elements()) {
        const int dofs[4] = {2 * e.node1(), 2 * e.node1() + 1, 2 * e.node2(), 2 * e.node2() + 1};
        const truss::Matrix4 k = e.stiffness();
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) K_ref[dofs[a]][dofs[b]] += k[a][b];
    }
    double max_diff = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            max_diff = std::max(max_diff, std::abs(model.stiffness()[i][j] - K_ref[i][j]));
    CHECK(max_diff < 1e-9);

    Model gauss = parse(text.str());
    const truss::SolveInfo gauss_info = gauss.solve(SolverType::Gauss);
    Model cg = parse(text.str());
    const truss::SolveInfo cg_info = cg.solve(SolverType::CG);
    CHECK(gauss_info.residual < 1e-10);
    CHECK(cg_info.residual < 1e-9);

    double u_max = 0.0, diff = 0.0;
    for (int d = 0; d < n; ++d) {
        u_max = std::max(u_max, std::abs(gauss.displacements()[d]));
        diff = std::max(diff, std::abs(gauss.displacements()[d] - cg.displacements()[d]));
    }
    CHECK(u_max > 0.0);
    CHECK(diff <= 1e-8 * u_max);

    double total_load_x = 0.0;
    for (int i = 0; i < nx; ++i) total_load_x += 1000.0 * i;
    CHECK_CLOSE(gauss.force_balance_x(), 0.0, 1e-6 * total_load_x);
    CHECK_CLOSE(gauss.force_balance_y(), 0.0, 1e-6 * 5000.0 * nx);
}

void test_mechanism_is_detected() {
    // A rectangle without a diagonal sways sideways: no unique solution.
    const std::string text = R"(
4 3 1 2
0 0 0
1 1000 0
2 0 1000
3 1000 1000
0 0 2 200000 100
1 2 3 200000 100
2 1 3 200000 100
2 1000 0
0 1 1
1 1 1
)";
    CHECK_THROWS_WITH(parse(text).solve(SolverType::Gauss), "singular at node");
    CHECK_THROWS_WITH(parse(text).solve(SolverType::CG), "mechanism");

    // No supports at all.
    CHECK_THROWS_WITH(parse("2 1 0 0\n0 0 0\n1 1 0\n0 0 1 1 1\n").solve(SolverType::Gauss),
                      "singular");
}

void test_comments_and_free_format() {
    Model model = parse(R"(
# a comment line
2 1   1 2   # counts
0 0 0       # node 0
1 1000      # values may continue
  0         # on the next line
0 0 1 1000 10
1 500 0
0 1 1
1 0 1
)");
    model.solve(SolverType::Gauss);
    CHECK_CLOSE(model.displacements()[2], 500.0 * 1000.0 / (1000.0 * 10.0), 1e-12);
}

void test_loads_and_supports_accumulate() {
    // Two loads on one node add up; a second support line cannot release a DOF.
    Model model = parse(
        "2 1 2 3\n0 0 0\n1 100 0\n0 0 1 1000 1\n1 3 0\n1 2 0\n0 1 1\n1 0 1\n1 0 0\n");
    CHECK(model.is_fixed(3));
    model.solve(SolverType::Gauss);
    CHECK_CLOSE(model.axial_forces()[0], 5.0, 1e-12);
    CHECK_CLOSE(model.reactions()[0], -5.0, 1e-12);
}

void test_input_errors() {
    const std::string tail = "0 0 1 1 1\n0 1 1\n1 0 1\n";  // element + supports

    CHECK_THROWS_WITH(Model::from_file("no_such_file.txt"), "cannot open");
    CHECK_THROWS_WITH(parse(""), "unexpected end of input");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n0 1 0\n" + tail), "test:3: node 0 is defined twice");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n5 1 0\n" + tail), "node 5 does not exist");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n1 1 abc\n" + tail), "expected a number");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n1 0 0\n" + tail), "test:4: element 0 has zero length");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n1 1 0\n0 0 1 -5 1\n0 1 1\n1 0 1\n"),
                      "E of element 0 must be positive");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n1 1 0\n0 0 0 1 1\n0 1 1\n1 0 1\n"),
                      "connects a node to itself");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n1 1 0\n0 0 1 1 1\n0 1 2\n1 0 1\n"),
                      "must be 0 (free) or 1 (fixed)");
    CHECK_THROWS_WITH(parse("2 1 0 2\n0 0 0\n1 1 0\n" + tail + "7\n"), "test:7: unexpected data");
    CHECK_THROWS_WITH(parse("1 1 0 0\n"), "at least 2 nodes");
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests = {
        {"single_bar", test_single_bar},
        {"two_bar_truss", test_two_bar_truss},
        {"four_node_truss_reference", test_four_node_truss_reference},
        {"element_stiffness", test_element_stiffness},
        {"assembly_symmetric_and_balanced", test_assembly_symmetric_and_balanced},
        {"large_grid_parallel_paths", test_large_grid_parallel_paths},
        {"mechanism_is_detected", test_mechanism_is_detected},
        {"comments_and_free_format", test_comments_and_free_format},
        {"loads_and_supports_accumulate", test_loads_and_supports_accumulate},
        {"input_errors", test_input_errors},
    };

    for (const auto& [name, test] : tests) {
        const int failures_before = g_failures;
        try {
            test();
        } catch (const std::exception& e) {
            report_failure(__FILE__, __LINE__, std::string(name) + " threw: " + e.what());
        }
        std::cout << (g_failures == failures_before ? "[ OK ] " : "[FAIL] ") << name << '\n';
    }
    std::cout << '\n' << g_checks << " checks, " << g_failures << " failures\n";
    return g_failures == 0 ? 0 : 1;
}
