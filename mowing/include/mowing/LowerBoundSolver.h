#ifndef TSPN_LOWERBOUNDSSOLVER_H
#define TSPN_LOWERBOUNDSSOLVER_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <random>
#include <utility>
#include <chrono>
#include <CGAL/ch_graham_andrew.h>

#include "utils/utils.hpp"
#include "utils/tours.hpp"
#include "utils/polygons.hpp"
#include "definitions.h"
#include "close_enough_tsp/CETSPSolver.h"
#include "ExactOffsetCalculator.h"
#include "mowing_utils.h"
#include "TSPSolver.h"
#include "tour_utils.h"


namespace mowing {
    class LowerBoundSolver {
    public:
        static const int INITIAL_STRATEGY_CH = 1;
        static const int INITIAL_STRATEGY_VERTICES = 2;

        static const int FOLLOWUP_STRATEGY_RANDOM = 3;
        static const int FOLLOWUP_STRATEGY_GRID = 4;
        static const int FOLLOWUP_STRATEGY_SKELETON = 5;

        typedef std::vector<Point> PointVector;
        typedef ExactOffsetCalculator::Gps_traits::General_polygon_2 ConicPolygon;
        typedef ExactOffsetCalculator::ConicPolygonVector ConicPolygonVector;

        struct lower_bound_solution {
            double cetsp_time;
            double tsp_time;
            int strategy;
            PointVector witnesses;
            PointVector k_means_centroids;
            PointVector base_witnesses;
            LowerBoundSolver::PointVector tour;
            LowerBoundSolver::PointVector tour_after_tsp;
            bool optimal;
            double lower_bound;
        };

        struct solution {
            Polygon_2 polygon;
            double radius;
            double lower_bound;
            std::vector<lower_bound_solution> lb_solutions;
        };

        LowerBoundSolver(Polygon_2 &polygon,
                         int initial_strategy, int followup_strategy,
                         double radius, double time,
                         std::size_t max_witness_size_initial,
                         std::size_t max_witness_size,
                         std::size_t max_iterations);

        solution solve();

        double getLowerBound();

    protected:
        Polygon_2 straight_line_polygon;

        std::unique_ptr<mowing::ExactOffsetCalculator> offset_calculator;

        double radius;
        int initial_strategy;
        int followup_strategy;

        Kernel::FT lower_bound;

        double time;

        std::size_t max_witness_size;
        std::size_t max_witness_size_initial;
        std::size_t max_iterations;

        void initializeOffsetCalculator();

        ConicPolygonVector computeUncoveredRegions(PointVector &tour);

        void getOptimalCETSPTour(lower_bound_solution &lb_solution);

        void addGridWitnesses(std::vector<Point> &witnesses, const ConicPolygon &P);
        void addRandomWitnesses(std::vector<Point> &witnesses, const ConicPolygon &P, std::size_t count);

        void addToWitnessSet(std::vector<Point> &witnesses, const ConicPolygon &P);

        void addInteriorWitnesses(std::vector<Point> &witnesses, CGAL::Polygon_2<Epick> &P, bool add_edges);

        void cleanUpAndAddWitnesses(std::vector<Point> &points, std::vector<Point> &witnesses);

        Polygon_2 approximatePolygon(const ConicPolygon &P);

        void updateLowerBound(double lb);
        virtual void updateLowerBound(const Kernel::FT &lb);

        void calculateDispersionApproximation(std::vector<Point> &witnesses, std::size_t n, Polygon_2* P = nullptr);
        PointVector kMeansSelection(std::vector<Point> &witnesses, std::size_t n, ConicPolygon &P);
    };
}

#endif
