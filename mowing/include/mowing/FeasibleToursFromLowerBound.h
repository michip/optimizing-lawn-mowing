#ifndef FEASIBLE_TOURS_FROM_LB_H
#define FEASIBLE_TOURS_FROM_LB_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <utility>
#include <chrono>

#include "utils/utils.hpp"
#include "utils/tours.hpp"
#include "definitions.h"
#include "close_enough_tsp/CETSPSolver.h"
#include "ExactOffsetCalculator.h"
#include "mowing_utils.h"
#include "TSPSolver.h"
#include "tour_utils.h"
#include "mowing/LowerBoundSolver.h"

namespace mowing {
    class FeasibleToursFromLowerBound {
    public:
        typedef std::vector<Point> PointVector;
        typedef ExactOffsetCalculator::ConicPolygonVector ConicPolygonVector;
        typedef ExactOffsetCalculator::Polygon_with_holes_2 Conic_Polygon_with_holes_2;
        typedef ExactOffsetCalculator::Gps_traits::General_polygon_2 ConicPolygon;
        typedef LowerBoundSolver::solution LBSolution;

        struct solution_step {
            FeasibleToursFromLowerBound::PointVector witness_set;
            FeasibleToursFromLowerBound::PointVector tour;
            std::vector<Conic_Polygon_with_holes_2> uncovered_area;
            double time;
        };

        struct solution_iteration {
            std::vector<solution_step> steps;
            double lower_bound;
            double upper_bound;
            double time;
        };

        struct solution {
            std::vector<solution_iteration> iterations;

            Polygon_2 polygon;
            double radius;
            double lower_bound;
            double upper_bound;
            bool solved_optimally;
            double time;
        };

        FeasibleToursFromLowerBound(LBSolution &solution, double time, std::size_t max_witness_size);

        void registerCallback(std::function<void(solution &)> &cb);

        solution solve();
        void fixSolution(solution &result);

        double getLowerBound();

        double getUpperBound();


    protected:
        Polygon_2 straight_line_polygon;
        LBSolution lower_bound_solution;

        std::unique_ptr<mowing::ExactOffsetCalculator> offset_calculator;
        std::function<void(solution &)> callback;

        double radius;

        Kernel::FT lower_bound;
        Kernel::FT upper_bound;

        double time;

        std::size_t max_witness_size;

        void initializeOffsetCalculator();

        virtual ConicPolygonVector computeUncoveredRegions(PointVector &tour);

        virtual void
        addTourForRegion(PointVector &tour, Conic_Polygon_with_holes_2 &region, PointVector &all_witnesses);

        virtual PointVector computeTourForConic(const ConicPolygon &P,
                                                PointVector &all_witnesses, Point *start_point);

        CETSPSolver::solution getOptimalCETSPTour(PointVector &witnesses, Point *start_point = nullptr);

        void addToWitnessSet(std::vector<Point> &witnesses, const ConicPolygon &P);

        void addInteriorWitnesses(std::vector<Point> &witnesses, CGAL::Polygon_2<Epick> &P, bool add_edges);

        void cleanUpAndAddWitnesses(std::vector<Point> &points, std::vector<Point> &witnesses);

        Polygon_2 approximatePolygon(const ConicPolygon &P);

        void updateLowerBound(double lb);
        virtual void updateLowerBound(const Kernel::FT &lb);
        void updateLowerBound(PointVector &solution);
        virtual void updateUpperBound(PointVector &solution);

        Kernel::FT computeGap();

        PointVector kMeansSelection(std::vector<Point> &witnesses, std::size_t n, const ConicPolygon &P);

        void writeBoundsAndExport(solution &sol, solution_iteration &iteration);
    };
}

#endif
