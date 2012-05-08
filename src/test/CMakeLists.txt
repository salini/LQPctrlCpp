

add_executable( test_minimal
    test_minimal.cpp
)
target_link_libraries( test_minimal
    Solver
)


add_executable( test_readFromFile
    test_readFromFile.cpp
)
target_link_libraries( test_readFromFile
    Solver
)


add_executable( test_readFromArray
    test_readFromArray.cpp
)
target_link_libraries( test_readFromArray
    Solver
)

#add_executable( test_random
#    test_random.cpp
#)
#target_link_libraries( test_random
#    Solver
#)
