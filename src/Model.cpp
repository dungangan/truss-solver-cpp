#include "truss/Model.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace truss {

namespace {

// ---------------------------------------------------------------------------
// Input reading
// ---------------------------------------------------------------------------

// Splits the input into whitespace separated tokens and remembers the line of each
// one, so errors can point at the right place. '#' starts a comment.
class TokenReader {
public:
    TokenReader(std::istream& in, std::string source) : source_{std::move(source)} {
        std::string line;
        int line_number = 0;
        while (std::getline(in, line)) {
            ++line_number;
            const std::size_t hash = line.find('#');
            if (hash != std::string::npos) {
                line.erase(hash);
            }
            std::istringstream words(line);
            std::string word;
            while (words >> word) {
                tokens_.push_back({word, line_number});
            }
        }
    }

    int next_int(const std::string& what) {
        const Token& token = next(what);
        const char* begin = token.text.c_str();
        char* end = nullptr;
        errno = 0;
        const long value = std::strtol(begin, &end, 10);
        if (end == begin || *end != '\0' || errno == ERANGE || value < INT_MIN ||
            value > INT_MAX) {
            fail("expected an integer for " + what + ", got '" + token.text + "'");
        }
        return static_cast<int>(value);
    }

    double next_double(const std::string& what) {
        const Token& token = next(what);
        const char* begin = token.text.c_str();
        char* end = nullptr;
        errno = 0;
        const double value = std::strtod(begin, &end);
        if (end == begin || *end != '\0' || errno == ERANGE || !std::isfinite(value)) {
            fail("expected a number for " + what + ", got '" + token.text + "'");
        }
        return value;
    }

    // Error at the most recently read token.
    [[noreturn]] void fail(const std::string& message) const {
        fail_at(pos_ > 0 ? tokens_[pos_ - 1].line : 1, message);
    }

    [[noreturn]] void fail_at(int line, const std::string& message) const {
        throw std::runtime_error(source_ + ":" + std::to_string(line) + ": " + message);
    }

    void expect_end() const {
        if (pos_ < tokens_.size()) {
            fail_at(tokens_[pos_].line, "unexpected data '" + tokens_[pos_].text +
                                            "' after the last support (check the counts "
                                            "on the first line)");
        }
    }

    int line() const { return pos_ > 0 ? tokens_[pos_ - 1].line : 1; }

private:
    struct Token {
        std::string text;
        int line;
    };

    const Token& next(const std::string& what) {
        if (pos_ >= tokens_.size()) {
            throw std::runtime_error(source_ + ": unexpected end of input, expected " + what);
        }
        return tokens_[pos_++];
    }

    std::string source_;
    std::vector<Token> tokens_;
    std::size_t pos_ = 0;
};

std::string str(int value) { return std::to_string(value); }

// Restores the formatting flags of a stream when it goes out of scope.
class StreamStateGuard {
public:
    explicit StreamStateGuard(std::ostream& os)
        : os_{os}, flags_{os.flags()}, precision_{os.precision()} {}
    ~StreamStateGuard() {
        os_.flags(flags_);
        os_.precision(precision_);
    }
    StreamStateGuard(const StreamStateGuard&) = delete;
    StreamStateGuard& operator=(const StreamStateGuard&) = delete;

private:
    std::ostream& os_;
    std::ios_base::fmtflags flags_;
    std::streamsize precision_;
};

constexpr int kMaxPrintedDofs = 30;  // larger matrices are not printed

// Adding +0.0 turns -0.0 into 0.0, so zero entries never print as "-0.0000e+00".
double printable(double value) { return value + 0.0; }

}  // namespace

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

Model Model::from_file(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open model file '" + path + "'");
    }
    return from_stream(file, path);
}

Model Model::from_stream(std::istream& in, const std::string& source_name) {
    TokenReader reader(in, source_name);
    Model model;

    const int n_nodes = reader.next_int("the number of nodes");
    if (n_nodes < 2) reader.fail("a truss needs at least 2 nodes");
    const int n_elements = reader.next_int("the number of elements");
    if (n_elements < 1) reader.fail("a truss needs at least 1 element");
    const int n_loads = reader.next_int("the number of loads");
    if (n_loads < 0) reader.fail("the number of loads cannot be negative");
    const int n_supports = reader.next_int("the number of supports");
    if (n_supports < 0) reader.fail("the number of supports cannot be negative");

    auto check_node = [&](int node) {
        if (node < 0 || node >= n_nodes) {
            reader.fail("node " + str(node) + " does not exist (valid ids: 0.." +
                        str(n_nodes - 1) + ")");
        }
    };

    // Nodes: "id x y". Any order, but every id 0..n-1 exactly once.
    std::vector<double> xs(n_nodes), ys(n_nodes);
    std::vector<char> node_seen(n_nodes, 0);
    for (int i = 0; i < n_nodes; ++i) {
        const int id = reader.next_int("a node id");
        check_node(id);
        if (node_seen[id]) reader.fail("node " + str(id) + " is defined twice");
        node_seen[id] = 1;
        xs[id] = reader.next_double("the x coordinate of node " + str(id));
        ys[id] = reader.next_double("the y coordinate of node " + str(id));
    }
    model.nodes_.reserve(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        model.nodes_.emplace_back(i, xs[i], ys[i]);
    }

    // Elements: "id node1 node2 E A".
    struct ElementInput {
        int node1, node2;
        double E, A;
        int line;
    };
    std::vector<ElementInput> element_input(n_elements);
    std::vector<char> element_seen(n_elements, 0);
    for (int i = 0; i < n_elements; ++i) {
        const int id = reader.next_int("an element id");
        if (id < 0 || id >= n_elements) {
            reader.fail("element id " + str(id) + " is out of range 0.." + str(n_elements - 1));
        }
        if (element_seen[id]) reader.fail("element " + str(id) + " is defined twice");
        element_seen[id] = 1;

        ElementInput& e = element_input[id];
        e.line = reader.line();
        e.node1 = reader.next_int("the first node of element " + str(id));
        check_node(e.node1);
        e.node2 = reader.next_int("the second node of element " + str(id));
        check_node(e.node2);
        if (e.node1 == e.node2) reader.fail("element " + str(id) + " connects a node to itself");
        e.E = reader.next_double("E of element " + str(id));
        if (!(e.E > 0.0)) reader.fail("E of element " + str(id) + " must be positive");
        e.A = reader.next_double("A of element " + str(id));
        if (!(e.A > 0.0)) reader.fail("A of element " + str(id) + " must be positive");
    }
    model.elements_.reserve(n_elements);
    for (int id = 0; id < n_elements; ++id) {
        const ElementInput& e = element_input[id];
        try {
            model.elements_.emplace_back(id, model.nodes_[e.node1], model.nodes_[e.node2], e.E,
                                         e.A);
        } catch (const std::invalid_argument& error) {
            reader.fail_at(e.line, error.what());
        }
        model.nodes_[e.node1].add_element(id);
        model.nodes_[e.node2].add_element(id);
    }

    // Loads: "node Fx Fy". Several loads on one node add up.
    model.f_.assign(2 * n_nodes, 0.0);
    for (int i = 0; i < n_loads; ++i) {
        const int node = reader.next_int("the node of a load");
        check_node(node);
        model.f_[2 * node] += reader.next_double("Fx of the load on node " + str(node));
        model.f_[2 * node + 1] += reader.next_double("Fy of the load on node " + str(node));
    }

    // Supports: "node fix_x fix_y" with 1 = fixed, 0 = free.
    model.fixed_.assign(2 * n_nodes, 0);
    for (int i = 0; i < n_supports; ++i) {
        const int node = reader.next_int("the node of a support");
        check_node(node);
        for (int direction = 0; direction < 2; ++direction) {
            const std::string name = direction == 0 ? "fix_x" : "fix_y";
            const int flag = reader.next_int(name + " of the support on node " + str(node));
            if (flag != 0 && flag != 1) reader.fail(name + " must be 0 (free) or 1 (fixed)");
            if (flag == 1) model.fixed_[2 * node + direction] = 1;
        }
    }

    reader.expect_end();
    return model;
}

int Model::num_free_dofs() const {
    int free = 0;
    for (const char fixed : fixed_) {
        if (!fixed) ++free;
    }
    return free;
}

// ---------------------------------------------------------------------------
// Assembly and solution
// ---------------------------------------------------------------------------

void Model::assemble() {
    const int n = num_nodes();
    const int ndof = num_dofs();
    K_.assign(ndof, std::vector<double>(ndof, 0.0));

    // Row-wise assembly. The iteration for node i writes only rows 2i and 2i+1:
    // for every element i-j it adds k to block (i,i) and -k to block (i,j).
    // Each row is written by exactly one thread, so no locks or atomics are needed.
#pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i) {
        std::vector<double>& row_x = K_[2 * i];
        std::vector<double>& row_y = K_[2 * i + 1];
        for (const int e : nodes_[i].elements()) {
            const T2D2& element = elements_[e];
            const Block2 k = element.stiffness_block();
            const int j = element.node1() == i ? element.node2() : element.node1();
            for (int b = 0; b < 2; ++b) {
                row_x[2 * i + b] += k[0][b];
                row_x[2 * j + b] -= k[0][b];
                row_y[2 * i + b] += k[1][b];
                row_y[2 * j + b] -= k[1][b];
            }
        }
    }
    solved_ = false;
}

SolveInfo Model::solve(SolverType type) {
    if (K_.empty()) {
        assemble();
    }
    const int ndof = num_dofs();

    // Supports: zero the row and column of every fixed DOF and keep its diagonal
    // entry. This forces u = 0 there, keeps the matrix symmetric positive definite
    // and does not change its scaling.
    Matrix K = K_;
    std::vector<double> f = f_;
    for (int d = 0; d < ndof; ++d) {
        if (!fixed_[d]) continue;
        const double diagonal = K[d][d] > 0.0 ? K[d][d] : 1.0;
        for (int j = 0; j < ndof; ++j) {
            K[d][j] = 0.0;
            K[j][d] = 0.0;
        }
        K[d][d] = diagonal;
        f[d] = 0.0;
    }

    SolveInfo info;
    info.type = type;
    try {
        if (type == SolverType::Gauss) {
            u_ = gauss_solve(std::move(K), f);
        } else {
            u_ = cg_solve(K, f, 1e-10, 10 * ndof + 100, &info.iterations);
        }
    } catch (const SingularMatrixError& error) {
        const int node = error.row() / 2;
        const char* direction = error.row() % 2 == 0 ? "x" : "y";
        throw std::runtime_error("the stiffness matrix is singular at node " + str(node) + " (" +
                                 direction +
                                 " direction): the truss is a mechanism or is not "
                                 "supported enough");
    }

    // Residual of the free equations. u is exactly zero at the fixed DOFs, so the
    // assembled K_ gives the same residual as the modified matrix.
    const std::vector<double> Ku = multiply(K_, u_);
    double r2 = 0.0;
    double f2 = 0.0;
    for (int d = 0; d < ndof; ++d) {
        if (fixed_[d]) continue;
        r2 += (Ku[d] - f[d]) * (Ku[d] - f[d]);
        f2 += f[d] * f[d];
    }
    info.residual = f2 > 0.0 ? std::sqrt(r2 / f2) : std::sqrt(r2);

    compute_results();
    solved_ = true;
    return info;
}

void Model::compute_results() {
    axial_forces_.assign(num_elements(), 0.0);

    // Sum the element end forces at every DOF; this equals K u. Element e pulls on
    // node1 with -N (c, s) and on node2 with +N (c, s).
    std::vector<double> internal(num_dofs(), 0.0);
    for (const T2D2& element : elements_) {
        const double N = element.axial_force(u_);
        axial_forces_[element.id()] = N;
        const int a = element.node1();
        const int b = element.node2();
        internal[2 * a] -= N * element.cos();
        internal[2 * a + 1] -= N * element.sin();
        internal[2 * b] += N * element.cos();
        internal[2 * b + 1] += N * element.sin();
    }

    // At a support, the part of K u that the applied load does not balance is the
    // reaction force.
    reactions_.assign(num_dofs(), 0.0);
    for (int d = 0; d < num_dofs(); ++d) {
        if (fixed_[d]) reactions_[d] = internal[d] - f_[d];
    }
}

double Model::force_balance_x() const {
    double sum = 0.0;
    for (int i = 0; i < num_nodes(); ++i) {
        sum += f_[2 * i] + (solved_ ? reactions_[2 * i] : 0.0);
    }
    return sum;
}

double Model::force_balance_y() const {
    double sum = 0.0;
    for (int i = 0; i < num_nodes(); ++i) {
        sum += f_[2 * i + 1] + (solved_ ? reactions_[2 * i + 1] : 0.0);
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

void Model::print_input(std::ostream& os) const {
    StreamStateGuard guard(os);
    os << std::fixed << std::setprecision(4);

    os << "Nodes\n"
       << "  node             x             y  support             Fx             Fy\n";
    for (const Node& node : nodes_) {
        const int i = node.id();
        const std::string support =
            std::string(fixed_[2 * i] ? "x" : "-") + (fixed_[2 * i + 1] ? "y" : "-");
        os << std::setw(6) << i << std::setw(14) << node.x() << std::setw(14) << node.y()
           << std::setw(9) << support << std::setw(15) << f_[2 * i] << std::setw(15)
           << f_[2 * i + 1] << '\n';
    }

    os << "\nElements\n"
       << "  elem  nodes             E             A        length\n";
    for (const T2D2& e : elements_) {
        const std::string nodes = str(e.node1()) + "-" + str(e.node2());
        os << std::setw(6) << e.id() << std::setw(7) << nodes << std::setw(14) << e.E()
           << std::setw(14) << e.A() << std::setw(14) << e.length() << '\n';
    }
}

void Model::print_results(std::ostream& os) const {
    if (!solved_) {
        os << "(model is not solved yet)\n";
        return;
    }
    StreamStateGuard guard(os);
    os << std::scientific << std::setprecision(6);

    os << "Nodal displacements\n"
       << "  node            ux            uy\n";
    for (int i = 0; i < num_nodes(); ++i) {
        os << std::setw(6) << i << std::setw(14) << printable(u_[2 * i]) << std::setw(14)
           << printable(u_[2 * i + 1]) << '\n';
    }

    double max_force = 0.0;
    for (const double N : axial_forces_) max_force = std::max(max_force, std::abs(N));

    os << "\nElement results (N > 0: tension)\n"
       << "  elem  nodes   axial force N    stress N/A  state\n";
    for (const T2D2& e : elements_) {
        const double N = axial_forces_[e.id()];
        const char* state = std::abs(N) <= 1e-9 * max_force ? "zero"
                            : N > 0.0                      ? "tension"
                                                           : "compression";
        const std::string nodes = str(e.node1()) + "-" + str(e.node2());
        os << std::setw(6) << e.id() << std::setw(7) << nodes << std::setw(16) << printable(N)
           << std::setw(14) << printable(N / e.A()) << "  " << state << '\n';
    }

    os << "\nSupport reactions\n"
       << "  node            Rx            Ry\n";
    for (int i = 0; i < num_nodes(); ++i) {
        if (!fixed_[2 * i] && !fixed_[2 * i + 1]) continue;
        os << std::setw(6) << i;
        for (int direction = 0; direction < 2; ++direction) {
            const int d = 2 * i + direction;
            if (fixed_[d]) {
                os << std::setw(14) << printable(reactions_[d]);
            } else {
                os << std::setw(14) << "-";
            }
        }
        os << '\n';
    }

    os << std::setprecision(3) << "\nForce balance (loads + reactions): sum Fx = "
       << force_balance_x() << ", sum Fy = " << force_balance_y() << '\n';
}

void Model::print_matrices(std::ostream& os) const {
    StreamStateGuard guard(os);
    os << std::scientific << std::setprecision(4);

    if (num_elements() <= kMaxPrintedDofs) {
        for (const T2D2& e : elements_) {
            os << "Element " << e.id() << " stiffness (DOFs " << 2 * e.node1() << ", "
               << 2 * e.node1() + 1 << ", " << 2 * e.node2() << ", " << 2 * e.node2() + 1
               << ")\n";
            for (const auto& row : e.stiffness()) {
                for (const double value : row) os << std::setw(13) << printable(value);
                os << '\n';
            }
        }
    } else {
        os << "Element stiffness matrices: skipped (more than " << kMaxPrintedDofs
           << " elements)\n";
    }

    os << '\n';
    if (K_.empty()) {
        os << "Global stiffness: not assembled\n";
    } else if (num_dofs() <= kMaxPrintedDofs) {
        os << "Global stiffness K (before supports)\n";
        for (const auto& row : K_) {
            for (const double value : row) os << std::setw(13) << printable(value);
            os << '\n';
        }
    } else {
        os << "Global stiffness: skipped (" << num_dofs() << " x " << num_dofs()
           << " is too large to print)\n";
    }
}

}  // namespace truss
