

add_library( abstractSolver
    abstractSolver.cpp
)

target_link_libraries( abstractSolver
    Solver_CGAL
    Solver_QuadProg
    Solver_COPL_QP
)


################################################################################
#
# Create CGAL Library
#
################################################################################
add_library( Solver_CGAL
    Solver_CGAL.cpp
)
target_link_libraries( Solver_CGAL
    abstractSolver
    ${CGAL_LIB}
    ${GMPXX_LIBRARIES}
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


#################################################################################
##
## Create COPL_QP Library
##
#################################################################################
add_library( Solver_COPL_QP
    Solver_COPL_QP.cpp
)
target_link_libraries( Solver_COPL_QP
    abstractSolver
    COPL_QP
)


#################################################################################
##
## Create OOQP Library
##
#################################################################################
#add_library( Solver_OOQP
#    Solver_OOQP.cpp
#)
#target_link_libraries( Solver_OOQP
#)

