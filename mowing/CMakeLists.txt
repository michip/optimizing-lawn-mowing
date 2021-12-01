add_library(mowing
        src/mowing_utils.cpp
        src/tour_utils.cpp
        src/ExactOffsetCalculator.cpp
        src/TSPSolver.cpp
        src/LowerBoundSolver.cpp
        include/mowing/Christofides.h
        src/Christofides.cpp
        src/FeketeApproximation.cpp
        src/FeasibleToursFromLowerBound.cpp
        include/mowing/tour_utils.h
        src/tour_utils.cpp)
set_target_properties(mowing PROPERTIES LINKER_LANGUAGE CXX)
util_setup_target(mowing
        LIBRARIES close_enough_tsp ${GUROBI} ${_GUROBI_CPP_STATIC_LIB} ${_GUROBI_DLL}
        INCLUDE_DIRS ${INCLUDE_DIRS}
        ${CGAL_INCLUDE_DIR}
        ${_GUROBI_CPP_INCLUDE_DIR}
        ${CMAKE_BASE_DIR}/include
        ${CMAKE_BASE_DIR}/utils/include
        ${CMAKE_CURRENT_LIST_DIR}/include)
set_target_properties(mowing PROPERTIES LINKER_LANGUAGE CXX)