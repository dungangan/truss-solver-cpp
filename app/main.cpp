// truss_solver: command line front end for the 2D truss finite element solver.

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "truss/Model.h"

#ifndef TRUSS_VERSION
#define TRUSS_VERSION "dev"
#endif

namespace {

struct Options {
    std::string model_path;
    truss::SolverType solver = truss::SolverType::Gauss;
    int threads = 0;  // 0: OpenMP default (OMP_NUM_THREADS, otherwise all cores)
    bool verbose = false;
    bool help = false;
};

class UsageError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

void print_usage(std::ostream& os) {
    os << "Usage: truss_solver <model.txt> [options]\n"
          "\n"
          "Options:\n"
          "  --solver gauss|cg   linear solver (default: gauss)\n"
          "  --threads N         number of OpenMP threads (default: OMP_NUM_THREADS\n"
          "                      or all cores)\n"
          "  -v, --verbose       also print the input tables and stiffness matrices\n"
          "  -h, --help          show this help\n";
}

Options parse_args(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value_of = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) throw UsageError(name + " needs a value");
            return argv[++i];
        };

        if (arg == "-h" || arg == "--help") {
            options.help = true;
        } else if (arg == "-v" || arg == "--verbose") {
            options.verbose = true;
        } else if (arg == "--solver") {
            const std::string name = value_of(arg);
            if (name == "gauss") {
                options.solver = truss::SolverType::Gauss;
            } else if (name == "cg") {
                options.solver = truss::SolverType::CG;
            } else {
                throw UsageError("unknown solver '" + name + "' (use gauss or cg)");
            }
        } else if (arg == "--threads") {
            const std::string text = value_of(arg);
            char* end = nullptr;
            const long n = std::strtol(text.c_str(), &end, 10);
            if (end == text.c_str() || *end != '\0' || n < 1 || n > 1024) {
                throw UsageError("--threads needs a positive integer, got '" + text + "'");
            }
            options.threads = static_cast<int>(n);
        } else if (!arg.empty() && arg[0] == '-') {
            throw UsageError("unknown option '" + arg + "'");
        } else if (options.model_path.empty()) {
            options.model_path = arg;
        } else {
            throw UsageError("only one model file can be given");
        }
    }
    if (!options.help && options.model_path.empty()) {
        throw UsageError("no model file given");
    }
    return options;
}

using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parse_args(argc, argv);
    } catch (const UsageError& error) {
        std::cerr << "error: " << error.what() << "\n\n";
        print_usage(std::cerr);
        return 2;
    }
    if (options.help) {
        print_usage(std::cout);
        return 0;
    }

    try {
#ifdef _OPENMP
        if (options.threads > 0) omp_set_num_threads(options.threads);
        const std::string threads = std::to_string(omp_get_max_threads());
#else
        if (options.threads > 1) {
            std::cerr << "warning: built without OpenMP, --threads is ignored\n";
        }
        const std::string threads = "1 (built without OpenMP)";
#endif

        Clock::time_point start = Clock::now();
        truss::Model model = truss::Model::from_file(options.model_path);
        const double t_read = ms_since(start);

        start = Clock::now();
        model.assemble();
        const double t_assemble = ms_since(start);

        start = Clock::now();
        const truss::SolveInfo info = model.solve(options.solver);
        const double t_solve = ms_since(start);

        std::cout << "truss_solver " << TRUSS_VERSION << '\n'
                  << "Model:   " << options.model_path << '\n'
                  << "Size:    " << model.num_nodes() << " nodes, " << model.num_elements()
                  << " elements, " << model.num_dofs() << " DOFs (" << model.num_free_dofs()
                  << " free)\n"
                  << "Solver:  " << truss::to_string(info.type);
        if (info.type == truss::SolverType::CG) {
            std::cout << ", " << info.iterations << " iterations";
        }
        std::cout << ", relative residual " << std::scientific << std::setprecision(2)
                  << info.residual << std::defaultfloat << '\n'
                  << "Threads: " << threads << "\n\n";

        if (options.verbose) {
            model.print_input(std::cout);
            std::cout << '\n';
            model.print_matrices(std::cout);
            std::cout << '\n';
        }
        model.print_results(std::cout);

        std::cout << std::fixed << std::setprecision(3) << "\nTiming [ms]: read " << t_read
                  << ", assemble " << t_assemble << ", solve " << t_solve << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
