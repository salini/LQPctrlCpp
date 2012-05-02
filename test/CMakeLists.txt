


add_executable( test_readFromFile
    test_readFromFile.cpp
)
target_link_libraries( test_readFromFile 
    genericSolver
)

add_executable( test_Solvers
    test_Solvers.cpp
)
target_link_libraries( test_Solvers
    genericSolver
)
