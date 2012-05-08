
################################################################################
#
# Generate solver libraries or dummy libraries if solver cannot be found
#
################################################################################
add_subdirectory(Solvers)
add_subdirectory(test)

add_library( Solver
    Solver.cpp
)

target_link_libraries( Solver
    abstractSolver
)

add_library( Constraint
    Constraint.cpp
)

add_library( Task
    Task.cpp
)