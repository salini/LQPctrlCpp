################################################################################
#
# Configuration
#
################################################################################
cmake_minimum_required(VERSION 2.8)
project(LQPctrlCpp)

set(CMAKE_MODULE_PATH      "${LQPctrlCpp_SOURCE_DIR}/cmake")

set(EXECUTABLE_OUTPUT_PATH "${LQPctrlCpp_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH    "${LQPctrlCpp_BINARY_DIR}/lib")


################################################################################
#
# Search external lib
#
################################################################################
##### Search for CGAL package
Find_Package(CGAL)
if(CGAL_FOUND)
    message(STATUS "CGAL package found.")
    link_directories(${CGAL_LIBRARIES_DIR})
    include_directories(${CGAL_INCLUDE_DIRS})
    add_definitions(-DCGAL_IS_AVAILABLE)
else()
    message(STATUS "The CGAL package has not been found.")
endif()

Find_Package(GMP)
if(GMP_FOUND)
    message(STATUS "GMP package found.")
    add_definitions(-DUSE_GMP_WIH_CGAL)
else()
    message(STATUS "The GMP package has not been found.")
endif()

##### The QuadProg package is available in the extra folder
add_definitions(-DQUADPROG_IS_AVAILABLE)


##### Search for Mosek package
if(MOSEK_FOUND)
    message(STATUS "MOSEK package found.")
    add_definitions(-DMOSEK_IS_AVAILABLE)
else()
    message(STATUS "The Mosek package has not been found.")
endif()


################################################################################
#
# Set include folders and create lib and executable
#
################################################################################

# Set include folder
include_directories(include)
include_directories(extra/quadprog-1.1)

add_subdirectory(test)
add_subdirectory(extra)
add_subdirectory(src)



