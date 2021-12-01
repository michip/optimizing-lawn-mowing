cmake_minimum_required(VERSION 3.16)
project(lawn_mowing)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_OSX_ARCHITECTURES x86_64)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # not necessary, but encouraged
set(CMAKE_C_FLAGS "-frounding-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -frounding-math")
set(target_link_options "-lgurobi_c++ -lgurobi60, --as-needed -arch;x86_64")

set(UTIL_BOOST_LIBRARIES system filesystem thread graph program_options)

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/util.cmake" NO_POLICY_SCOPE)

set(_CPLEX_CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})

set(CMAKE_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

find_package(GUROBI REQUIRED)
find_package(CPLEX REQUIRED)
if (APPLE)
    set(CPLEX_INCLUDE_DIRS /usr/local/include/)
elseif(UNIX)
    set(CPLEX_INCLUDE_DIRS "/opt/tmp/ibm/ILOG/CPLEX_Studio129/cplex/include/ /opt/tmp/ibm/ILOG/CPLEX_Studio129/concert/include/")
endif ()
find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(JsonCpp)

add_subdirectory(close_enough_tsp)
add_subdirectory(mowing)

message(STATUS "CGAL found at ${CGAL_LIBRARY}")
message(STATUS "CGAL include dir at ${CGAL_INCLUDE_DIR}")

add_executable(lower_bounds lower_bounds.cpp)
util_setup_target(lower_bounds
        LIBRARIES
        mowing ${JsonCpp_LIBRARIES}
        INCLUDE_DIRS
        ${CPLEX_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/close_enough_tsp/include
        ${CMAKE_CURRENT_SOURCE_DIR}/mowing/include
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/include
        ${JsonCpp_INCLUDE_DIRS}
        ${CGAL_INCLUDE_DIR})
set_target_properties(lower_bounds PROPERTIES LINKER_LANGUAGE CXX)


add_executable(approximation approximation.cpp)
util_setup_target(approximation
        LIBRARIES
        mowing ${JsonCpp_LIBRARIES}
        INCLUDE_DIRS
        ${CPLEX_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/close_enough_tsp/include
        ${CMAKE_CURRENT_SOURCE_DIR}/mowing/include
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/include
        ${JsonCpp_INCLUDE_DIRS}
        ${CGAL_INCLUDE_DIR})
set_target_properties(approximation PROPERTIES LINKER_LANGUAGE CXX)


add_executable(cetsp_solver_from_lb solve_from_lb_solution.cpp)
util_setup_target(cetsp_solver_from_lb
        LIBRARIES
        mowing ${JsonCpp_LIBRARIES}
        INCLUDE_DIRS
        ${CPLEX_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/close_enough_tsp/include
        ${CMAKE_CURRENT_SOURCE_DIR}/mowing/include
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/include
        ${JsonCpp_INCLUDE_DIRS}
        ${CGAL_INCLUDE_DIR})
set_target_properties(cetsp_solver_from_lb PROPERTIES LINKER_LANGUAGE CXX)


add_executable(fix_solution fix_solution.cpp)
util_setup_target(fix_solution
        LIBRARIES
        mowing ${JsonCpp_LIBRARIES}
        INCLUDE_DIRS
        ${CPLEX_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/close_enough_tsp/include
        ${CMAKE_CURRENT_SOURCE_DIR}/mowing/include
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/include
        ${JsonCpp_INCLUDE_DIRS}
        ${CGAL_INCLUDE_DIR})
set_target_properties(fix_solution PROPERTIES LINKER_LANGUAGE CXX)