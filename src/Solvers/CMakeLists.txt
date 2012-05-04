

add_library( abstractSolver
    abstractSolver.cpp
)

target_link_libraries( abstractSolver
    Solver_QuadProg
)


#################################################################################
##
## Create QuadProg++ Library
##
#################################################################################
add_library( Solver_QuadProg
    Solver_QuadProg.cpp
)
target_link_libraries( Solver_QuadProg
    abstractSolver
    QuadProg
)

#################################################################################
##
## Create Mosek Library
##
#################################################################################
#add_library( Solver_Mosek
#    Solver_Mosek.cpp
#)
#target_link_libraries( Solver_Mosek
#)




