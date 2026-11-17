cmake_minimum_required(VERSION 3.16)

project(truss_solver VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(TRUSS_USE_OPENMP "Parallelise assembly and solvers with OpenMP" ON)
option(TRUSS_BUILD_TESTS "Build the test suite" ON)

function(truss_set_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

# --- library -------------------------------------------------------------
add_library(truss STATIC
  src/T2D2.cpp
  src/Model.cpp
  src/Solvers.cpp
)
target_include_directories(truss PUBLIC ${PROJECT_SOURCE_DIR}/include)
truss_set_warnings(truss)

if(TRUSS_USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(truss PUBLIC OpenMP::OpenMP_CXX)
  else()
    message(WARNING "OpenMP not found: building the serial version")
  endif()
endif()
if(NOT TARGET OpenMP::OpenMP_CXX AND NOT MSVC)
  # Without OpenMP the '#pragma omp' lines are ignored; do not warn about them.
  target_compile_options(truss PRIVATE -Wno-unknown-pragmas)
endif()

# --- command line program ------------------------------------------------
add_executable(truss_solver app/main.cpp)
target_link_libraries(truss_solver PRIVATE truss)
target_compile_definitions(truss_solver PRIVATE TRUSS_VERSION="${PROJECT_VERSION}")
truss_set_warnings(truss_solver)

# --- tests ---------------------------------------------------------------
if(TRUSS_BUILD_TESTS)
  enable_testing()

  add_executable(truss_tests tests/test_truss.cpp)
  target_link_libraries(truss_tests PRIVATE truss)
  truss_set_warnings(truss_tests)
  add_test(NAME unit_tests COMMAND truss_tests)

  set(example ${PROJECT_SOURCE_DIR}/examples/four_node_truss.txt)
  add_test(NAME example_gauss COMMAND truss_solver ${example} --solver gauss)
  add_test(NAME example_cg COMMAND truss_solver ${example} --solver cg --threads 2)
  add_test(NAME cli_missing_file COMMAND truss_solver does_not_exist.txt)
  set_tests_properties(cli_missing_file PROPERTIES WILL_FAIL TRUE)
endif()
