#include "mowing/LowerBoundSolver.h"

namespace mowing {


    LowerBoundSolver::LowerBoundSolver(Polygon_2 &polygon, int initial_strategy, int followup_strategy,
                                       double radius, double time, std::size_t max_witness_size_initial,
                                       std::size_t max_witness_size,
                                       std::size_t max_iterations) {
        this->straight_line_polygon = polygon;

        if (this->straight_line_polygon.is_clockwise_oriented()) {
            this->straight_line_polygon.reverse_orientation();
        }

        this->radius = radius;
        this->initializeOffsetCalculator();

        this->lower_bound = 0;
        this->initial_strategy = initial_strategy;
        this->followup_strategy = followup_strategy;

        this->max_witness_size = max_witness_size;
        this->max_witness_size_initial = max_witness_size_initial;
        this->max_iterations = max_iterations;

        this->time = time;
    }

    void LowerBoundSolver::initializeOffsetCalculator() {
        this->offset_calculator = std::make_unique<mowing::ExactOffsetCalculator>(this->straight_line_polygon,
                                                                                  this->radius + 1e-3, false);
    }

    LowerBoundSolver::ConicPolygonVector
    LowerBoundSolver::computeUncoveredRegions(PointVector &tour) {
        this->offset_calculator->computeUncoveredRegions(tour);

        auto result = ExactOffsetCalculator::ConicPolygonVector();
        this->offset_calculator->polygon_set.polygons_with_holes(std::back_inserter(result));
        return result;
    }

    LowerBoundSolver::solution LowerBoundSolver::solve() {
        auto result = solution{this->straight_line_polygon,
                               this->radius,
                               CGAL::to_double(this->lower_bound),
                               std::vector<lower_bound_solution>()};

        try {

            auto all_witnesses = PointVector(); // Holds the witnesses for all iterations

            for (std::size_t i = 0; i < this->max_iterations; i++) {

                auto centroids = PointVector();
                auto base_witnesses = PointVector();

                int strategy = i == 0 ? this->initial_strategy : this->followup_strategy;

                if (i == 0) {
                    switch (this->initial_strategy) {
                        case INITIAL_STRATEGY_VERTICES:
                            this->addToWitnessSet(all_witnesses, this->offset_calculator->base_polygon);
                            break;
                        case INITIAL_STRATEGY_CH:
                            CGAL::ch_graham_andrew(this->straight_line_polygon.vertices_begin(),
                                                   this->straight_line_polygon.vertices_end(),
                                                   std::back_inserter(all_witnesses));
                            break;
                        default:
                            throw std::runtime_error(
                                    "The given initial_strategy " + std::to_string(this->initial_strategy) +
                                    " is not implemented");
                    }

                    base_witnesses.insert(base_witnesses.end(), all_witnesses.begin(), all_witnesses.end());
                    if (all_witnesses.size() > this->max_witness_size_initial) {
                        this->calculateDispersionApproximation(all_witnesses, this->max_witness_size_initial);
                    }
                }
                else {
                    auto tour = result.lb_solutions.back().tour_after_tsp.empty() ? result.lb_solutions.back().tour :
                                result.lb_solutions.back().tour_after_tsp;

                    auto uncovered = this->computeUncoveredRegions(tour);

                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::vector<double> probabilities;
                    std::map<int, int> counts;

                    int idx = 0;
                    for (auto &region: uncovered) {
                        counts[idx++] = 0;
                        Polygon_2 poly = this->approximatePolygon(region.outer_boundary());

                        if (poly.is_simple()) {
                            probabilities.emplace_back(CGAL::to_double(CGAL::abs(poly.area())));
                        } else if(!probabilities.empty()) {
                            probabilities.emplace_back(
                                    std::accumulate(probabilities.begin(), probabilities.end(), 0.0) /
                                            (double) probabilities.size());
                        } else {
                            probabilities.emplace_back(1);
                        }
                    }
                    std::discrete_distribution<> d(probabilities.begin(), probabilities.end());
                    for (std::size_t n = 0; n < this->max_witness_size; ++n) {
                        ++counts[d(gen)];
                    }

                    idx = 0;
                    for (auto &region: uncovered) {
                        auto witnesses = PointVector();

                        std::size_t count = counts[idx++];

                        auto boundary = region.outer_boundary();

                        switch (this->followup_strategy) {
                            case FOLLOWUP_STRATEGY_SKELETON:
                                this->addToWitnessSet(witnesses, boundary);
                                break;
                            case FOLLOWUP_STRATEGY_GRID:
                                this->addGridWitnesses(witnesses, boundary);
                                break;
                            case FOLLOWUP_STRATEGY_RANDOM:
                                this->addRandomWitnesses(witnesses, boundary, count);
                                break;
                            default:
                                throw std::runtime_error(
                                        "The given strategy " + std::to_string(this->followup_strategy) +
                                        " is not implemented");
                        }

                        base_witnesses.insert(base_witnesses.end(), witnesses.begin(), witnesses.end());

                        if (witnesses.size() > count) {
                            auto current_centroids = this->kMeansSelection(witnesses, count, boundary);
                            centroids.insert(centroids.end(), current_centroids.begin(), current_centroids.end());
                        }

                        all_witnesses.insert(all_witnesses.end(), witnesses.begin(), witnesses.end());
                    }

                    std::cout << all_witnesses.size() << " witness size" << std::endl;

                }

                result.lb_solutions.emplace_back(lower_bound_solution{
                        0, 0,
                        strategy,
                        PointVector(all_witnesses),
                        PointVector(centroids),
                        PointVector(base_witnesses),
                        PointVector(),
                        PointVector(), false, 0});
                this->getOptimalCETSPTour(result.lb_solutions.back());
                result.lower_bound = this->getLowerBound();
                this->offset_calculator->initializeUncoveredRegions();
            }
        } catch (std::runtime_error &ex) {
            result.lower_bound = this->getLowerBound();
            std::cout << "An error " << ex.what() << std::endl;
            return result;
        } catch (GRBException &ex) {
            result.lower_bound = this->getLowerBound();
            std::cout << "An error " << ex.getErrorCode() << ": " << ex.getMessage() << std::endl;
            return result;
        }

        return result;
    }

    void
    LowerBoundSolver::addInteriorWitnesses(std::vector<Point> &witnesses, CGAL::Polygon_2<Epick> &P, bool add_edges) {
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

    void LowerBoundSolver::addGridWitnesses(std::vector<Point> &witnesses, const ConicPolygon &P) {
        Polygon_2 polygon = this->approximatePolygon(P);
        std::vector<Point> points = std::vector<Point>();

        auto bbox = P.bbox();

        double minX = bbox.xmin(), maxX = bbox.xmax(),
                minY = bbox.ymin(), maxY = bbox.ymax();

        for (double x = minX; x <= maxX; x += this->radius * 0.25) {
            for (double y = minY; y <= maxY; y += this->radius * 0.25) {
                auto p = Point(x, y);
                if (polygon.has_on_bounded_side(p) == CGAL::ON_BOUNDED_SIDE) {
                    points.emplace_back(p);
                }
            }
        }

        this->cleanUpAndAddWitnesses(points, witnesses);
    }

    void LowerBoundSolver::addRandomWitnesses(std::vector<Point> &witnesses, const ConicPolygon &P, std::size_t count) {
        Polygon_2 polygon = this->approximatePolygon(P);

        for (std::size_t i = 0; i < count; i++) {
            witnesses.emplace_back(::utils::random_point_inside_polygon(polygon));
        }
    }

    void LowerBoundSolver::addToWitnessSet(std::vector<Point> &witnesses, const ConicPolygon &P) {
        std::vector<Point> points = std::vector<Point>();
        CGAL::Polygon_2<Epick> poly;
        Polygon_2 exactPolygon = this->approximatePolygon(P);
        for (auto it = exactPolygon.vertices_begin(); it != exactPolygon.vertices_end(); it++) {
            poly.push_back(Epick::Point_2(CGAL::to_double(it->x()), CGAL::to_double(it->y().approx())));
            if (&this->offset_calculator->base_polygon == &P) points.emplace_back(*it);
        }

        std::size_t invalid_count = 0;

        auto add_point = [&poly, &points, &invalid_count, this](Epick::Point_2 p) {
            bool valid = true;
            for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++) {
                if (*it == p) {
                    valid = false;
                    invalid_count++;
                    break;
                }
            }

            if (valid) points.emplace_back(p.x(), p.y());
        };

        if (&this->offset_calculator->base_polygon != &P) {
            if (poly.is_simple()) {
                //this->addInteriorWitnesses(points, poly, true);

                auto skeleton = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

                for (auto it = skeleton->halfedges_begin(); it != skeleton->halfedges_end(); it++) {
                    if (it->is_bisector()) {
                        add_point(Epick::Point_2(it->opposite()->vertex()->point().x(),
                                                 it->opposite()->vertex()->point().y()));
                        add_point(Epick::Point_2(it->vertex()->point().x(), it->vertex()->point().y()));
                    }
                }

                mowing::utils::get_points_along_half_edges_fixed_size(skeleton->halfedges_begin(),
                                                                      skeleton->halfedges_end(),
                                                                      5, points);
            } else {
                mowing::utils::get_points_along_edges(poly.edges_begin(),
                                                      poly.edges_end(), this->radius / 5, points);
            }
        }

        this->cleanUpAndAddWitnesses(points, witnesses);
    }

    void LowerBoundSolver::getOptimalCETSPTour(lower_bound_solution &lb_solution) {

        using std::chrono::high_resolution_clock;
        using std::chrono::duration;
        using std::chrono::milliseconds;

        auto t1 = high_resolution_clock::now();
        auto solver = CETSPSolver(lb_solution.witnesses, this->radius, this->time);

        auto solution = solver.solve();
        auto t2 = high_resolution_clock::now();

        lb_solution.cetsp_time = ((duration<double, std::milli>) (t2 - t1)).count();
        lb_solution.optimal = solution.optimal_solution_found;

        auto &tour = solution.points;

        if (tour.empty()) {
            this->updateLowerBound(solution.lower_bound);
            lb_solution.lower_bound = solution.lower_bound;
            throw std::runtime_error("Empty tour generation");
        }

        auto tour_end = tour.back();
        mowing::utils::cleanup_tour(tour);

        if (tour.empty()) { // All tour points were the same
            tour.emplace_back(tour_end);
        }

        lb_solution.tour = tour;

        if (solution.optimal_solution_found) {
            auto length = ::utils::compute_tour_length(tour);
            this->updateLowerBound(length);
            lb_solution.lower_bound = CGAL::to_double(length);
        } else {
            this->updateLowerBound(solution.lower_bound);
            lb_solution.lower_bound = solution.lower_bound;
        }

        if (!solution.optimal_solution_found && tour.size() > 2) {
            t1 = high_resolution_clock::now();
            auto new_witnesses = PointVector(lb_solution.witnesses);
            auto tsp_solver = mowing::tsp::TSPSolver(new_witnesses, tour, this->time);
            auto new_tour = tsp_solver.solve().points;
            t2 = high_resolution_clock::now();

            std::cout << "Improved cleanup tour (maybe) old " << ::utils::compute_tour_length(tour) << " vs. new "
                      << ::utils::compute_tour_length(new_tour) << std::endl;

            if (::utils::compute_tour_length(tour) > ::utils::compute_tour_length(new_tour)) {
                lb_solution.tsp_time = ((duration<double, std::milli>) (t2 - t1)).count();
                lb_solution.tour_after_tsp = tour;
            }
        }

    }


    void LowerBoundSolver::updateLowerBound(const Kernel::FT &lb) {
        this->lower_bound = CGAL::max(lb, this->lower_bound);
    }

    void LowerBoundSolver::updateLowerBound(double lb) {
        this->updateLowerBound((Kernel::FT) lb);
    }

    double LowerBoundSolver::getLowerBound() {
        return CGAL::to_double(this->lower_bound);
    }

    void LowerBoundSolver::cleanUpAndAddWitnesses(std::vector<Point> &points, std::vector<Point> &witnesses) {
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

    Polygon_2 LowerBoundSolver::approximatePolygon(const ConicPolygon &P) {
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
        } while(!polygon.is_simple() && iteration++ < maxIterations);

        return polygon;
    }

    void LowerBoundSolver::calculateDispersionApproximation(std::vector<Point> &witnesses, std::size_t n,
                                                            Polygon_2 *P) {
        std::vector<double> distances_to_boundary;
        for (auto &witness : witnesses) {
            if (P == nullptr) distances_to_boundary.emplace_back(0);
            else {
                std::vector<double> distances;
                for (auto it = P->edges_begin(); it != P->edges_end(); it++) {
                    auto s = mowing::utils::minimum_distance(it->source(), it->target(), witness);
                    distances.emplace_back(CGAL::to_double(s->squared_length()));
                }
                distances_to_boundary.emplace_back(*std::min_element(distances.begin(), distances.end()));
            }
        }

        using IndexPair = std::pair<std::size_t, std::size_t>;
        std::map<IndexPair, double> distances;

        std::set<size_t> selected_indices;
        std::vector<size_t> unselected_indices;

        double maxDistance = 0;
        double maxDistanceToBoundary = *std::max_element(distances_to_boundary.begin(), distances_to_boundary.end());

        for (std::size_t i = 0; i < witnesses.size(); i++) {
            unselected_indices.emplace_back(i);
            for (std::size_t j = i + 1; j < witnesses.size(); j++) {
                distances[std::make_pair(i, j)] = CGAL::to_double(CGAL::squared_distance(witnesses[i],
                                                                                         witnesses[j]));
                maxDistance = std::max(maxDistance, distances[std::make_pair(i, j)]);
            }
        }

        if (P != nullptr) {
            for (std::size_t i = 0; i < witnesses.size(); i++) {
                for (std::size_t j = i + 1; j < witnesses.size(); j++) {
                    distances[std::make_pair(i, j)] = distances[std::make_pair(i, j)] / maxDistance +
                                                      +0.5 * (distances_to_boundary[i] / maxDistanceToBoundary +
                                                              distances_to_boundary[j] / maxDistanceToBoundary);
                }
            }
        }
        // Do 2-Approximation of p-dispersion
        auto longestDistancePair = std::max_element(distances.begin(), distances.end(),
                                                    [&distances_to_boundary](const std::pair<IndexPair, double> &p1,
                                                                             const std::pair<IndexPair, double> &p2) {
                                                        return p1.second < p2.second;
                                                    });

        auto delete_from_unselected = [&unselected_indices](const std::size_t &idx) {
            auto position = std::find(unselected_indices.begin(), unselected_indices.end(), idx);
            if (position != unselected_indices.end())
                unselected_indices.erase(position);
        };

        selected_indices.insert(longestDistancePair->first.first);
        selected_indices.insert(longestDistancePair->first.second);

        delete_from_unselected(longestDistancePair->first.first);
        delete_from_unselected(longestDistancePair->first.second);


        while (selected_indices.size() < n) {
            double overall_max_min_distance = 0;
            auto index_to_add = unselected_indices.front();

            for (auto &i: unselected_indices) {
                double min_distance = longestDistancePair->second;
                for (auto &j: selected_indices) {
                    min_distance = std::min(min_distance, distances[std::make_pair(std::min(i, j), std::max(i, j))]);
                }

                if (overall_max_min_distance < min_distance) {
                    index_to_add = i;
                    overall_max_min_distance = min_distance;
                }
            }

            selected_indices.insert(index_to_add);
            delete_from_unselected(index_to_add);
        }

        std::sort(unselected_indices.begin(), unselected_indices.end());  // Make sure the container is sorted
        for (auto it = unselected_indices.rbegin(); it != unselected_indices.rend(); it++) {
            witnesses.erase(witnesses.begin() + (long) *it);
        }

    }

    LowerBoundSolver::PointVector LowerBoundSolver::kMeansSelection(std::vector<Point> &witnesses,
                                                                    std::size_t n, ConicPolygon &P) {
        std::vector<Point> centroids;

        if(n == 0) {
            witnesses.clear();
            return centroids;
        }

        std::vector<std::size_t> labels(witnesses.size(), 0);
        this->addRandomWitnesses(centroids, P, n); // Adding n random centroids

        assert(centroids.size() == n);

        std::vector<Point> oldCentroids;
        std::size_t currentIteration = 0;

        auto shouldStop = [&currentIteration, &centroids, &oldCentroids] () {
            bool stop = true;

            if(currentIteration > 200) return true;

            for(std::size_t i = 0; i < centroids.size(); i++) {
                if(centroids[i] != oldCentroids[i]) {
                    stop = false;
                    break;
                }
            }

            return stop;
        };

        auto getPointsForCentroid = [&witnesses, &labels](std::size_t k) {
            std::vector<Point> pointsForCentroid;
            for(std::size_t i = 0; i < labels.size(); i++) {
                if(labels[i] == k) {
                    pointsForCentroid.emplace_back(witnesses[i]);
                }
            }
            return pointsForCentroid;
        };

        auto calculateLabels = [&centroids, &witnesses, &labels] () {
            for(std::size_t i = 0; i < witnesses.size(); i++) {
                auto minDistance = CGAL::squared_distance(witnesses[i], centroids[0]);
                std::size_t label = 0;

                for(std::size_t k = 1; k < centroids.size(); k++) {
                    auto currentDistance = CGAL::squared_distance(witnesses[i], centroids[k]);
                    if(currentDistance < minDistance) {
                        minDistance = currentDistance;
                        label = k;
                    }
                }

                labels[i] = label;
            }
        };

        auto generateRandomPoint = [&P, this](){
            std::vector<Point> centroidDummy; // Doing this to reuse the function
            this->addRandomWitnesses(centroidDummy, P, 1); // Adding 1 random centroids
            return *centroidDummy.begin();
        };

        auto calculateCentroids = [&centroids, &getPointsForCentroid, &generateRandomPoint] () {
            for(std::size_t k = 0; k < centroids.size(); k++) {
                std::vector<Point> pointsForCentroid = getPointsForCentroid(k);
                if(pointsForCentroid.empty()) {
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
        } while(!shouldStop());

        std::vector<Point> newWitnessSet;
        std::default_random_engine re(std::random_device{}());
        for(std::size_t k = 0; k < centroids.size(); k++) {
            std::vector<Point> pointsForCentroid = getPointsForCentroid(k);
            // std::cout << "CENTROID " << k << " with " << pointsForCentroid.size() << " points" << std::endl;

            if(pointsForCentroid.empty()) {
                newWitnessSet.emplace_back(generateRandomPoint());
            } else {
                newWitnessSet.emplace_back(*::utils::select_randomly(pointsForCentroid.begin(), pointsForCentroid.end(), re));
            }
        }

        witnesses.clear();
        witnesses.insert(witnesses.end(), newWitnessSet.begin(), newWitnessSet.end());

        return centroids;
    }

}

