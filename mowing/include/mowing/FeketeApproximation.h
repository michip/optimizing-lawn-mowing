#ifndef TSPN_FEKETE_APPROXIMATION_H
#define TSPN_FEKETE_APPROXIMATION_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <utility>

#include "utils/utils.hpp"
#include "utils/tours.hpp"
#include "definitions.h"
#include "close_enough_tsp/CETSPSolver.h"
#include "ExactOffsetCalculator.h"
#include "mowing_utils.h"
#include "TSPSolver.h"
#include "tour_utils.h"
#include <CGAL/Boolean_set_operations_2.h>

namespace mowing {
    class FeketeApproximation {
    public:
        typedef std::vector<Point> PointVector;

        struct solution {
            Polygon_2 polygon;
            FeketeApproximation::PointVector points;
            FeketeApproximation::PointVector tour;
            double radius;
            double lower_bound;
            double upper_bound;
            bool solved_optimally;
            double runtime;
        };

        FeketeApproximation(Polygon_2 &polygon, double radius, double time);

        solution solve();

    protected:
        Polygon_2 straight_line_polygon;

        double radius;
        double time;

        std::vector<Point> generateHexagonalGrid();
    };
}

#endif
