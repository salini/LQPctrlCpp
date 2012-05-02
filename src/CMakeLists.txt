
################################################################################
#
# Generate solver libraries or dummy libraries if solver cannot be found
#
################################################################################
add_subdirectory(Solvers)

add_library( genericSolver
    genericSolver.cpp
)

target_link_libraries( genericSolver
    abstractSolver
)

