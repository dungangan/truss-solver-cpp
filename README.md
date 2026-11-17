# truss-solver-cpp

[![CI](https://github.com/dungangan/truss-solver-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/dungangan/truss-solver-cpp/actions/workflows/ci.yml)

A small finite element solver for 2D pin-jointed trusses, written in C++17 with
OpenMP. It reads a plain text model, assembles the global stiffness matrix in
parallel, solves for the nodal displacements and reports member forces, stresses
and support reactions.

The code is meant to be read: each step of the method has its own short function,
and the test suite checks the results against closed-form solutions.

## Features

- **T2D2 element**: 2-node linear truss element (named after the Abaqus element)
- **Parallel assembly** with OpenMP, written so that no two threads ever write to
  the same matrix entry (no locks or atomics)
- **Two linear solvers**: Gaussian elimination and the conjugate gradient method
- **Results**: displacements, axial forces, stresses, tension/compression,
  support reactions, and a force balance check
- **Clear errors**: input mistakes are reported with file name and line, and an
  unstable truss (a mechanism) is detected and named instead of giving garbage
- **Tests**: analytic cases, an independent reference solution, parallel vs
  serial assembly, input validation

## Build

Requirements: CMake 3.16+ and a C++17 compiler. OpenMP is used when available
(GCC, Clang with libomp, MSVC); otherwise a serial version is built.

```bash
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release
```

## Usage

```
truss_solver <model.txt> [options]

  --solver gauss|cg   linear solver (default: gauss)
  --threads N         number of OpenMP threads (default: OMP_NUM_THREADS or all cores)
  -v, --verbose       also print the input tables and stiffness matrices
  -h, --help          show this help
```

Example: [examples/four_node_truss.txt](examples/four_node_truss.txt), a braced
frame pinned at the bottom with a 50 kN horizontal load at the top left. It has
five members and four free DOFs, so it is statically indeterminate.

```
$ ./build/truss_solver examples/four_node_truss.txt
truss_solver 1.0.0
Model:   examples/four_node_truss.txt
Size:    4 nodes, 5 elements, 8 DOFs (4 free)
Solver:  gauss, relative residual 1.46e-16
Threads: 8

Nodal displacements
  node            ux            uy
     0  0.000000e+00  0.000000e+00
     1  0.000000e+00  0.000000e+00
     2  5.393521e+00  1.314208e+00
     3  4.866502e+00 -1.185792e+00

Element results (N > 0: tension)
  elem  nodes   axial force N    stress N/A  state
     0    0-2    3.942623e+04  1.752277e+02  tension
     1    2-3   -2.371584e+04 -1.054038e+02  compression
     2    0-3    4.275435e+04  1.900193e+02  tension
     3    2-1   -4.738444e+04 -2.105975e+02  compression
     4    1-3   -3.557377e+04 -1.581056e+02  compression

Support reactions
  node            Rx            Ry
     0 -2.371584e+04 -7.500000e+04
     1 -2.628416e+04  7.500000e+04

Force balance (loads + reactions): sum Fx = -5.821e-11, sum Fy = 1.455e-11

Timing [ms]: read 0.251, assemble 1.073, solve 0.109
```

## Input format

Values are separated by spaces or newlines; everything after `#` is a comment.
Use any consistent set of units (the example uses N, mm and MPa).

```
<nodes> <elements> <loads> <supports>

<node id> <x> <y>                     one line per node
<element id> <node1> <node2> <E> <A>  one line per element
<node> <Fx> <Fy>                      one line per load
<node> <fix_x> <fix_y>                one line per support, 1 = fixed, 0 = free
```

- Node ids are `0 .. nodes-1` and element ids are `0 .. elements-1`. They can be
  listed in any order, but each id must appear exactly once.
- Several loads on the same node add up.
- A pin is `1 1`, a roller is `0 1` or `1 0`. Prescribed displacements are zero.

## How it works

**Element stiffness.** For an element of length $L$ with direction
$c = \cos\theta$, $s = \sin\theta$:

$$
\mathbf{k} = \frac{EA}{L}\begin{bmatrix} c^2 & cs \\ cs & s^2 \end{bmatrix},
\qquad
\mathbf{K}_e = \begin{bmatrix} \mathbf{k} & -\mathbf{k} \\ -\mathbf{k} & \mathbf{k} \end{bmatrix}
$$

Length and direction are computed once when the element is created, so elements
are read-only during assembly.

**Parallel assembly.** Node $i$ owns the global DOFs $2i$ and $2i+1$. The loop
runs over nodes in parallel. For every element $i$–$j$ at node $i$ it adds
$\mathbf{k}$ to block $(i,i)$ and $-\mathbf{k}$ to block $(i,j)$ of $\mathbf{K}$.
The thread for node $i$ only writes rows $2i$ and $2i+1$, so the threads never
write to the same entry and need no synchronisation.

**Supports.** For each fixed DOF the row and column of $\mathbf{K}$ are set to
zero, the diagonal entry is kept and the load is set to zero. This forces
$u = 0$ there and keeps $\mathbf{K}$ symmetric positive definite.

**Solvers.**
- `gauss`: Gaussian elimination without pivoting, which is safe for a symmetric
  positive definite matrix. The rows below the pivot are updated in parallel, and
  rows that are already zero are skipped (truss matrices are banded). A zero pivot
  means the truss is a mechanism. The error names the node and direction.
- `cg`: conjugate gradient. Matrix-vector products and dot products run in
  parallel. It stops at relative residual $10^{-10}$.

**Results.** The axial force is $N = \frac{EA}{L}\left(c\,\Delta u_x + s\,\Delta u_y\right)$
(positive in tension) and the stress is $\sigma = N/A$. Adding up the element end
forces gives $\mathbf{K}\mathbf{u}$. At a support, the part not balanced by the
applied load is the reaction. Loads and reactions must add up to zero, and the
program prints this sum as a check.

## Verification

`truss_tests` (run by `ctest`) checks:

| Test | Reference |
|---|---|
| Single bar under tension | $u = FL/EA$, $N = F$ |
| Symmetric two-bar truss | $N = P / 2\sin\theta$, $\delta = PL / (2EA\sin^2\theta)$, reactions |
| Four-node example (indeterminate) | independent NumPy solution, both solvers, 10 significant digits |
| Element stiffness | symmetry, zero force under rigid translation, $EA/L$ under unit stretch |
| Assembled $\mathbf{K}$ | symmetric, rows sum to zero |
| 150-node braced grid | parallel assembly equals element-by-element assembly; Gauss and CG agree |
| Mechanisms | a frame without a diagonal and an unsupported bar are rejected |
| Input errors | duplicate or missing ids, zero length, bad numbers, wrong counts |

## Project layout

```
include/truss/   Node.h, T2D2.h, Model.h, Solvers.h
src/             T2D2.cpp, Model.cpp (input, assembly, supports, results), Solvers.cpp
app/main.cpp     command line program
tests/           test suite (no external framework)
examples/        example model
```

## Limitations

- Dense matrix storage: memory grows as $(2n)^2$, which suits models up to a few
  thousand nodes.
- Linear elastic material and small displacements. No buckling check for members
  in compression.
- Only nodal loads and zero prescribed displacements.
