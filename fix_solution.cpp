#include <stdlib.h>
#include <iostream>
#include <string>
#include <json/writer.h>
#include <json/json.h>
#include <algorithm>
#include "utils/utils.hpp"
#include "mowing/FeasibleToursFromLowerBound.h"
#include "mowing/LowerBoundSolver.h"
#include "mowing/ExactOffsetCalculator.h"
#include "utils/export.hpp"
#include "utils/import.hpp"
#include "utils/utils.hpp"



int main(int argc, char *argv[]) {

    if (argc != 7) {
        std::cout << "Too few or many arguments" << std::endl;
        return 0;
    }


    Polygon_2 polygon;

    std::string line;
    std::string file_name = argv[1];
    std::string out_file = argv[2];
    std::string in_file = argv[3];
    std::string in_file_solution = argv[4];
    double time = atof(argv[5]);
    auto max_witness_size = (std::size_t) stoul(argv[6]);

    std::ifstream input_file(file_name);
    if (!input_file.is_open()) {
        std::cout << "Could not open file " << file_name << std::endl;
        return 0;
    }
    std::getline(input_file, line);
    input_file.close();

    std::stringstream ss;
    ss << line;
    ss >> polygon;
    std::cout << "Input polygon " << file_name << " has " << polygon.size() << " inside_points" << std::endl;
    std::cout << "Input polygon is simple " << polygon.is_simple() << " and has area " << polygon.area() << std::endl;

    auto lb_solution = utils::importJsonLB(in_file);
    auto solution = utils::importSolutionJson(in_file_solution);

    std::cout << "Solution input " << lb_solution.polygon.size() << " has area " << lb_solution.polygon.area()
              << std::endl;

    if (CGAL::abs(lb_solution.polygon.area() - polygon.area()) > 0.01) {
        throw std::runtime_error("Polygons do not match");
    }

    try {
        auto solver = mowing::FeasibleToursFromLowerBound(lb_solution, time, max_witness_size);

        std::function<void(mowing::FeasibleToursFromLowerBound::solution &)> callback =
                [&out_file](mowing::FeasibleToursFromLowerBound::solution &sol) {
                    utils::exportJson(out_file, sol);
                };
        solver.registerCallback(callback);

        solver.fixSolution(solution);
        utils::exportJson(out_file, solution);

    } catch (std::exception &ex) {
        std::cout << ex.what() << std::endl;
    }
}
