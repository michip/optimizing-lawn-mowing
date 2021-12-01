#include "mowing/FeasibleToursFromLowerBound.h"


namespace mowing {


    FeasibleToursFromLowerBound::FeasibleToursFromLowerBound(LBSolution &solution, double time,
                                                             std::size_t max_witness_size) {
        this->straight_line_polygon = solution.polygon;

        if (this->straight_line_polygon.is_clockwise_oriented()) {
            this->straight_line_polygon.reverse_orientation();
        }

        this->lower_bound_solution = solution;

        this->radius = solution.radius;
        this->initializeOffsetCalculator();

        this->lower_bound = 0;
        this->upper_bound = 999999;

        this->max_witness_size = max_witness_size;

        this->time = time;
    }

    void FeasibleToursFromLowerBound::initializeOffsetCalculator() {
        this->offset_calculator = std::make_unique<mowing::ExactOffsetCalculator>(this->straight_line_polygon,
                                                                                  this->radius + 1e-3, false);
    }

    void FeasibleToursFromLowerBound::fixSolution(solution &result) {
        using std::chrono::high_resolution_clock;
        using std::chrono::duration;
        using std::chrono::milliseconds;


        bool solved;

        PointVector tour;
        auto all_witnesses = PointVector();

        auto full_start_time = high_resolution_clock::now();

        for (std::size_t i = 0; i < result.iterations.size(); i++) {
            auto &iteration = result.iterations[i];
            auto lb_solution = this->lower_bound_solution.lb_solutions[i];
            if (iteration.steps.empty()) continue;
            try {
                auto iteration_start_time = high_resolution_clock::now();
                this->updateLowerBound(iteration.lower_bound);

                tour = iteration.steps.back().tour;

                bool iterationSolved = this->computeUncoveredRegions(tour).empty();
                this->offset_calculator->initializeUncoveredRegions();

                if (iterationSolved) {
                    std::cout << "Iteration already solved" << std::endl;
                    continue;
                }

                result.time -= iteration.steps.back().time;
                iteration.time -= iteration.steps.back().time;
                iteration.steps.pop_back(); // Remove the last step so that it can be replaced.

                if (iteration.steps.empty()) { // Restart the step if it was not the first one. Else start over.
                    tour = lb_solution.tour;
                    all_witnesses = lb_solution.witnesses;
                } else {
                    all_witnesses = iteration.steps.back().witness_set;
                    tour = iteration.steps.back().tour;
                }


                do {
                    auto uncovered = this->computeUncoveredRegions(tour);
                    iteration.steps.push_back(solution_step{PointVector(all_witnesses),
                                                            PointVector(tour),
                                                            uncovered, 0});

                    auto step_start_time = high_resolution_clock::now();
                    solved = uncovered.empty();

                    if (solved) {
                        this->updateUpperBound(tour);
                    } else {
                        for (auto &region: uncovered) {
                            this->addTourForRegion(tour, region, all_witnesses);
                        }
                    }

                    auto step_end_time = high_resolution_clock::now();
                    iteration.steps.back().time = ((duration<double, std::milli>) (step_end_time -
                                                                                   step_start_time)).count();

                    std::cout << "Finished step for iteration " << i << std::endl;
                    this->writeBoundsAndExport(result, iteration);
                } while (!solved);

                std::cout << "Found feasible solution for iteration " << i << " opt between: "
                          << this->getLowerBound() << " and "
                          << this->getUpperBound() << " Gap: " <<
                          CGAL::to_double(this->computeGap()) <<
                          " solved opt " << (this->computeGap() < 1e-5)
                          << std::endl;

                this->offset_calculator->initializeUncoveredRegions();
                result.solved_optimally = this->computeGap() < 1e-5;

                auto iteration_end_time = high_resolution_clock::now();
                iteration.time += ((duration<double, std::milli>) (iteration_end_time - iteration_start_time)).count();

                this->writeBoundsAndExport(result, iteration);

            } catch (std::runtime_error &ex) {
                std::cout << "An error " << ex.what() << std::endl;
            }
        }

        result.lower_bound = this->getLowerBound();
        result.upper_bound = this->getUpperBound();

        auto full_end_time = high_resolution_clock::now();
        result.time += ((duration<double, std::milli>) (full_end_time - full_start_time)).count();
    }


    FeasibleToursFromLowerBound::solution FeasibleToursFromLowerBound::solve() {
        using std::chrono::high_resolution_clock;
        using std::chrono::duration;
        using std::chrono::milliseconds;

        auto result = solution{std::vector<solution_iteration>(),
                               this->straight_line_polygon,
                               this->radius,
                               CGAL::to_double(this->lower_bound),
                               CGAL::to_double(this->upper_bound),
                               false, 0};

        bool solved;

        PointVector tour;
        auto all_witnesses = PointVector();

        auto full_start_time = high_resolution_clock::now();

        for (auto &lb_solution: this->lower_bound_solution.lb_solutions) {
            try {
                tour = lb_solution.tour;
                all_witnesses = lb_solution.witnesses;

                result.iterations.push_back(solution_iteration{
                        std::vector<solution_step>(),
                        lb_solution.lower_bound,
                        999999, 0});

                auto iteration_start_time = high_resolution_clock::now();

                auto &iteration = result.iterations.back();
                this->updateLowerBound(lb_solution.lower_bound);

                do {
                    auto uncovered = this->computeUncoveredRegions(tour);
                    iteration.steps.push_back(solution_step{PointVector(all_witnesses),
                                                            PointVector(tour),
                                                            uncovered, 0});

                    auto step_start_time = high_resolution_clock::now();
                    solved = uncovered.empty();

                    if (solved) {
                        this->updateUpperBound(tour);
                    } else {
                        for (auto &region: uncovered) {
                            this->addTourForRegion(tour, region, all_witnesses);
                        }
                    }

                    auto step_end_time = high_resolution_clock::now();
                    iteration.steps.back().time = ((duration<double, std::milli>) (step_end_time -
                                                                                   step_start_time)).count();

                    std::cout << "Finished step" << std::endl;
                    this->writeBoundsAndExport(result, iteration);
                } while (!solved);

                std::cout << "Found feasible solution opt between: "
                          << this->getLowerBound() << " and "
                          << this->getUpperBound() << " Gap: " <<
                          CGAL::to_double(this->computeGap()) <<
                          " solved opt " << (this->computeGap() < 1e-5)
                          << std::endl;

                this->offset_calculator->initializeUncoveredRegions();
                result.solved_optimally = this->computeGap() < 1e-5;

                auto iteration_end_time = high_resolution_clock::now();
                iteration.time = ((duration<double, std::milli>) (iteration_end_time - iteration_start_time)).count();

                this->writeBoundsAndExport(result, iteration);

            } catch (std::runtime_error &ex) {
                std::cout << "An error " << ex.what() << std::endl;
            }
        }

        result.lower_bound = this->getLowerBound();
        result.upper_bound = this->getUpperBound();

        auto full_end_time = high_resolution_clock::now();
        result.time = ((duration<double, std::milli>) (full_end_time - full_start_time)).count();


        return result;
    }

    void FeasibleToursFromLowerBound::registerCallback(std::function<void(solution & )> &cb) {
        this->callback = cb;
    }

    void FeasibleToursFromLowerBound::writeBoundsAndExport(solution &sol, solution_iteration &iteration) {
        iteration.lower_bound = this->getLowerBound();
        iteration.upper_bound = this->getUpperBound();

        sol.lower_bound = iteration.lower_bound;
        sol.upper_bound = iteration.upper_bound;

        if (this->callback != nullptr) this->callback(sol);
    }

    Kernel::FT FeasibleToursFromLowerBound::computeGap() {
        return CGAL::abs(this->upper_bound - this->lower_bound) / this->upper_bound;
    }

    void
    FeasibleToursFromLowerBound::addInteriorWitnesses(std::vector<Point> &witnesses, CGAL::Polygon_2<Epick> &P,
                                                      bool add_edges) {
        if (P.is_clockwise_oriented()) {
            P.reverse_orientation();
        }

        if (add_edges) mowing::utils::get_points_along_edges(P.edges_begin(), P.edges_end(), this->radius, witnesses);

        for (auto it = P.vertices_begin(); it != P.vertices_end(); it++) {
            witnesses.emplace_back(it->x(), it->y());
        }

        /*
        auto offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(this->radius, P);

        for (auto &offset_polygon: offset_polygons) {
            this->addInteriorWitnesses(witnesses, *offset_polygon, false);
        }*/
    }

    void FeasibleToursFromLowerBound::addToWitnessSet(std::vector<Point> &witnesses, const ConicPolygon &P) {
        std::vector<Point> points = std::vector<Point>();

        CGAL::Polygon_2<Epick> poly;
        auto exactPolygon = this->approximatePolygon(P);
        for (auto it = exactPolygon.vertices_begin(); it != exactPolygon.vertices_end(); it++) {
            poly.push_back(Epick::Point_2(CGAL::to_double(it->x()), CGAL::to_double(it->y())));
            points.emplace_back(*it);
        }

        std::cout << "Equal Eval " << (&this->offset_calculator->base_polygon == &P) << std::endl;

        if (&this->offset_calculator->base_polygon != &P) {
            if (poly.is_simple()) {
                this->addInteriorWitnesses(points, poly, true);

                auto skeleton = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

                for (auto it = skeleton->halfedges_begin(); it != skeleton->halfedges_end(); it++) {
                    if (it->is_bisector()) {
                        points.emplace_back(it->opposite()->vertex()->point().x(),
                                            it->opposite()->vertex()->point().y());
                        points.emplace_back(it->vertex()->point().x(), it->vertex()->point().y());
                    }
                }

                mowing::utils::get_points_along_half_edges(skeleton->halfedges_begin(),
                                                           skeleton->halfedges_end(),
                                                           this->radius, points);
            } else {
                mowing::utils::get_points_along_edges(poly.edges_begin(), poly.edges_end(), this->radius, points);
            }
        }

        this->cleanUpAndAddWitnesses(points, witnesses);
    }

    void FeasibleToursFromLowerBound::addTourForRegion(PointVector &tour,
                                                       Conic_Polygon_with_holes_2 &region,
                                                       PointVector &all_witnesses) {

        auto shortest_connecting_result = mowing::utils::shortest_connecting_segment(tour, region);

        auto shortest_segment = shortest_connecting_result.first;
        auto edge_to_replace = shortest_connecting_result.second;

        Point start_point = shortest_segment->source();

        auto subtour = this->computeTourForConic(region.outer_boundary(),
                                                 all_witnesses,
                                                 &start_point);

        long start_index = 0;

        // Find index of the splitting point in the subtour (might be slightly different due to double precision)
        for (long i = 0; i < (long) subtour.size(); i++) {
            if (CGAL::has_smaller_distance_to_point(shortest_segment->source(), subtour[i], subtour[start_index])) {
                start_index = i;
            }
        }

        // Modify vector so that start vector is at front
        std::rotate(subtour.begin(), subtour.begin() + start_index, subtour.end());
        subtour[0] = start_point;

        for (auto it = tour.begin(); it != tour.end(); it++) {

            auto next = it + 1;
            if (next == tour.end()) next = tour.begin();

            if (edge_to_replace == nullptr) {
                if (*it == shortest_segment->source()) {
                    tour.insert(it, subtour.begin(), subtour.end());
                    break;
                }
            } else if (*it == edge_to_replace->first && *next == edge_to_replace->second) {
                subtour.emplace_back(start_point);
                tour.insert(next, subtour.begin(), subtour.end());
                break;
            }
        }

        mowing::utils::cleanup_tour(tour);
    }

    FeasibleToursFromLowerBound::PointVector
    FeasibleToursFromLowerBound::computeTourForConic(const ConicPolygon &P, PointVector &all_witnesses,
                                                     Point *start_point) {
        auto witnesses = PointVector();
        this->addToWitnessSet(witnesses, P);

        if (witnesses.size() > this->max_witness_size) {
            this->kMeansSelection(witnesses, this->max_witness_size, P);
        }

        auto solution = this->getOptimalCETSPTour(witnesses, start_point);
        auto &tour = solution.points;

        if (solution.optimal_solution_found) {
            this->updateLowerBound(tour);
        } else {
            this->updateLowerBound(solution.lower_bound);
        }

        all_witnesses.insert(all_witnesses.end(), witnesses.begin(), witnesses.end());

        return tour;
    }

    CETSPSolver::solution FeasibleToursFromLowerBound::getOptimalCETSPTour(PointVector &witnesses, Point *start_point) {
        auto solver = CETSPSolver(witnesses, *start_point, this->radius, this->time);

        std::streambuf *old = std::cout.rdbuf(); // <-- save
        std::stringstream ss;
        //std::cout.rdbuf(ss.rdbuf());       // <-- redirect

        auto solution = solver.solve();
        auto &tour = solution.points;
        std::cout.rdbuf(old);              // <-- restore

        if (tour.empty()) {
            this->updateLowerBound(solution.lower_bound);
            std::cout << "Could not generate a nonempty tour. Using TSP instead." << std::endl;
            auto new_witnesses = PointVector(witnesses);
            auto tsp_solver = mowing::tsp::TSPSolver(new_witnesses, this->time);
            tour = tsp_solver.solve().points;
            return solution;
        }

        auto tour_end = tour.back();
        mowing::utils::cleanup_tour(tour);

        if (tour.empty()) { // All tour points were the same
            tour.emplace_back(tour_end);
            //this->updateLowerBound(solution.lower_bound);
            //throw std::runtime_error("Empty tour after cleanup");
        }

        if (!solution.optimal_solution_found && tour.size() > 2) {
            auto new_witnesses = PointVector(witnesses);
            if (start_point != nullptr) new_witnesses.emplace_back(*start_point);
            auto tsp_solver = mowing::tsp::TSPSolver(new_witnesses, tour, this->time);
            auto new_tour = tsp_solver.solve().points;

            std::cout << "Improved cleanup tour (maybe) old " << ::utils::compute_tour_length(tour) << " vs. new "
                      << ::utils::compute_tour_length(new_tour) << std::endl;

            if (::utils::compute_tour_length(tour) > ::utils::compute_tour_length(new_tour)) {
                tour = new_tour;
            }
        }

        return solution;
    }

    FeasibleToursFromLowerBound::ConicPolygonVector
    FeasibleToursFromLowerBound::computeUncoveredRegions(PointVector &tour) {
        this->offset_calculator->computeUncoveredRegions(tour);

        auto result = ExactOffsetCalculator::ConicPolygonVector();
        this->offset_calculator->polygon_set.polygons_with_holes(std::back_inserter(result));
        return result;
    }

    void FeasibleToursFromLowerBound::updateLowerBound(const Kernel::FT &lb) {
        this->lower_bound = CGAL::max(lb, this->lower_bound);
    }

    void FeasibleToursFromLowerBound::updateLowerBound(double lb) {
        this->updateLowerBound((Kernel::FT) lb);
    }

    void FeasibleToursFromLowerBound::updateLowerBound(PointVector &solution) {
        this->updateLowerBound(::utils::compute_tour_length(solution));
    }

    void FeasibleToursFromLowerBound::updateUpperBound(PointVector &solution) {
        this->upper_bound = CGAL::min(::utils::compute_tour_length(solution), this->upper_bound);
    }

    double FeasibleToursFromLowerBound::getLowerBound() {
        return CGAL::to_double(this->lower_bound);
    }

    double FeasibleToursFromLowerBound::getUpperBound() {
        return CGAL::to_double(this->upper_bound);
    }

    void
    FeasibleToursFromLowerBound::cleanUpAndAddWitnesses(std::vector<Point> &points, std::vector<Point> &witnesses) {
        for (auto it = points.begin(); it != points.end();) {
            bool pointAlreadyPresent = false;

            for (auto it2 = it + 1; it2 != points.end(); it2++) {
                if (*it == *it2) {
                    pointAlreadyPresent = true;
                    break;
                }
            }

            for (auto &witness : witnesses) {
                if (*it == witness) {
                    pointAlreadyPresent = true;
                    break;
                }
            }

            if (pointAlreadyPresent) {
                it = points.erase(it);
            } else if (this->straight_line_polygon.has_on_unbounded_side(*it)) {
                auto segment = mowing::utils::shortest_connecting_segment(this->straight_line_polygon, *it);
                (*it) = segment->target();
                it++;
            } else {
                it++;
            }
        }


        for (auto &p: points) {
            witnesses.push_back(p);
        }
    }

    Polygon_2 FeasibleToursFromLowerBound::approximatePolygon(const ConicPolygon &P) {
        std::size_t resolution = 0;

        Polygon_2 polygon;
        std::size_t maxIterations = 10;
        std::size_t iteration = 0;

        do {
            resolution += 10;

            for (auto it = P.curves_begin(); it != P.curves_end(); it++) {
                auto approximation = std::vector<std::pair<double, double>>();

                it->polyline_approximation(resolution, std::back_inserter(approximation));
                if (!it->is_directed_right()) {
                    std::reverse(approximation.begin(), approximation.end());
                }

                approximation.pop_back();

                for (auto &pair: approximation) {
                    polygon.push_back(Point(pair.first, pair.second));
                }
            }
        } while (!polygon.is_simple() && iteration++ < maxIterations);

        return polygon;
    }

    FeasibleToursFromLowerBound::PointVector FeasibleToursFromLowerBound::kMeansSelection(std::vector<Point> &witnesses,
                                                                                          std::size_t n,
                                                                                          const ConicPolygon &P) {
        auto addRandomWitnesses = [this](std::vector<Point> &W, const ConicPolygon &Pol, std::size_t count) {
            Polygon_2 polygon = this->approximatePolygon(Pol);

            for (std::size_t i = 0; i < count; i++) {
                W.emplace_back(::utils::random_point_inside_polygon(polygon));
            }
        };

        std::vector<Point> centroids;

        if (n == 0) {
            witnesses.clear();
            return centroids;
        }

        std::vector<std::size_t> labels(witnesses.size(), 0);
        addRandomWitnesses(centroids, P, n); // Adding n random centroids

        assert(centroids.size() == n);

        std::vector<Point> oldCentroids;
        std::size_t currentIteration = 0;

        auto shouldStop = [&currentIteration, &centroids, &oldCentroids]() {
            bool stop = true;

            if (currentIteration > 200) return true;

            for (std::size_t i = 0; i < centroids.size(); i++) {
                if (centroids[i] != oldCentroids[i]) {
                    stop = false;
                    break;
                }
            }

            return stop;
        };

        auto getPointsForCentroid = [&witnesses, &labels](std::size_t k) {
            std::vector<Point> pointsForCentroid;
            for (std::size_t i = 0; i < labels.size(); i++) {
                if (labels[i] == k) {
                    pointsForCentroid.emplace_back(witnesses[i]);
                }
            }
            return pointsForCentroid;
        };

        auto calculateLabels = [&centroids, &witnesses, &labels]() {
            for (std::size_t i = 0; i < witnesses.size(); i++) {
                auto minDistance = CGAL::squared_distance(witnesses[i], centroids[0]);
                std::size_t label = 0;

                for (std::size_t k = 1; k < centroids.size(); k++) {
                    auto currentDistance = CGAL::squared_distance(witnesses[i], centroids[k]);
                    if (currentDistance < minDistance) {
                        minDistance = currentDistance;
                        label = k;
                    }
                }

                labels[i] = label;
            }
        };

        auto generateRandomPoint = [&P, &addRandomWitnesses]() {
            std::vector<Point> centroidDummy; // Doing this to reuse the function
            addRandomWitnesses(centroidDummy, P, 1); // Adding 1 random centroids
            return *centroidDummy.begin();
        };

        auto calculateCentroids = [&centroids, &getPointsForCentroid, &generateRandomPoint]() {
            for (std::size_t k = 0; k < centroids.size(); k++) {
                std::vector<Point> pointsForCentroid = getPointsForCentroid(k);
                if (pointsForCentroid.empty()) {
                    centroids[k] = generateRandomPoint();
                } else {
                    centroids[k] = CGAL::centroid(pointsForCentroid.begin(), pointsForCentroid.end());
                }
            }
        };

        do {
            oldCentroids = std::vector<Point>(centroids);
            calculateLabels();
            calculateCentroids();
            currentIteration++;
        } while (!shouldStop());

        std::vector<Point> newWitnessSet;
        std::default_random_engine re(std::random_device{}());
        for (std::size_t k = 0; k < centroids.size(); k++) {
            std::vector<Point> pointsForCentroid = getPointsForCentroid(k);
            if (pointsForCentroid.empty()) {
                newWitnessSet.emplace_back(generateRandomPoint());
            } else {
                newWitnessSet.emplace_back(
                        *::utils::select_randomly(pointsForCentroid.begin(), pointsForCentroid.end(), re));
            }
        }

        witnesses.clear();
        witnesses.insert(witnesses.end(), newWitnessSet.begin(), newWitnessSet.end());

        return centroids;
    }
}
