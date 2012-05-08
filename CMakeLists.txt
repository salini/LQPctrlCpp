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

##### The QuadProg package is available in the extra folder
add_definitions(-DQUADPROG_IS_AVAILABLE)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


##### Search for Mosek package
#if(MOSEK_FOUND)
#    message(STATUS "MOSEK package found.")
#    add_definitions(-DMOSEK_IS_AVAILABLE)
#else()
#    message(STATUS "The Mosek package has not been found.")
#endif()


################################################################################
#
# Set include folders and create lib and executable
#
################################################################################

# Set include folder
include_directories(include)
include_directories(extra/quadprog-1.2)

add_subdirectory(extra)
add_subdirectory(src)



