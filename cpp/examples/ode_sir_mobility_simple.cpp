
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/model.h"
#include "ode_sir_mobility/parameters.h"
#include "ode_sir_mobility/regions.h"
#include "ode_sir_mobility/contact_location.h"
#include "memilio/io/io.h"

mio::IOResult<std::vector<std::vector<std::vector<int>>>> read_path_mobility(const std::string& filename)
{
    BOOST_OUTCOME_TRY(num_lines, mio::count_lines(filename));

    if (num_lines == 0) {
        std::vector<std::vector<std::vector<int>>> arr(0, std::vector<std::vector<int>>(0));
        return mio::success(arr);
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(mio::StatusCode::FileNotFound, filename);
    }

    std::vector<std::vector<std::vector<int>>> arr(std::sqrt(num_lines),
                                                   std::vector<std::vector<int>>(std::sqrt(num_lines)));

    try {
        std::string tp;
        while (getline(file, tp)) {
            auto line  = mio::split(tp, ' ');
            int indx_x = std::stoi(line[0]);
            int indx_y = std::stoi(line[1]);
            if (indx_x != indx_y) {
                auto path = std::accumulate(line.begin() + 2, line.end(), std::string(""));

                // string -> vector of integers
                std::vector<int> path_vec;

                // Remove the square brackets and \r
                path = path.substr(1, path.size() - 3);
                std::stringstream ss(path);
                std::string token;

                // get numbers and save them in path_vec
                while (std::getline(ss, token, ',')) {
                    path_vec.push_back(std::stoi(token));
                }

                // Sorted by end location
                for (int number : path_vec) {
                    if (number != indx_x && number != indx_y) {
                        arr[indx_x][indx_y].push_back(number);
                    }
                }
            }
        }
    }
    catch (std::runtime_error& ex) {
        return failure(mio::StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return mio::success(arr);
}

mio::IOResult<void> run(const std::string& filename, mio::osirmobility::Model& model)
{
    BOOST_OUTCOME_TRY(mobility_paths, read_path_mobility(filename));
    size_t n_regions = (size_t)model.parameters.get_num_regions();
    for (size_t i = 0; i < n_regions; i++) {
        for (size_t j = 0; j < n_regions; j++) {
            if (j == i) {
                continue;
            }
            std::sort(mobility_paths[i][j].begin(), mobility_paths[i][j].end());
            std::vector<int> intersection_int;
            std::vector<mio::osirmobility::Region> intersection_region(intersection_int.size(),
                                                                       mio::osirmobility::Region(0));
            for (size_t k = 0; k < n_regions; k++) {
                if (k == i || k == j) {
                    continue;
                }
                std::sort(mobility_paths[k][j].begin(), mobility_paths[k][j].end());
                std::set_intersection(mobility_paths[i][j].begin(), mobility_paths[i][j].end(),
                                      mobility_paths[k][j].begin(), mobility_paths[k][j].end(),
                                      std::back_inserter(intersection_int));

                if (intersection_int.begin() != intersection_int.end()) {
                    intersection_region.push_back(mio::osirmobility::Region(k));
                    intersection_int.pop_back();
                }
            }
            if (intersection_region.begin() != intersection_region.end()) {
                model.parameters.get<mio::osirmobility::PathIntersections>()[mio::Index(
                    mio::osirmobility::Region(i), mio::osirmobility::Region(j))] = intersection_region;
            }
        }
    }
    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 50.;
    double dt   = 1;

    size_t number_regions              = 4;
    size_t number_age_groups           = 1;
    size_t total_population_per_region = 10;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& filename = "";

    mio::osirmobility::Model model(number_regions, number_age_groups);

    for (size_t i = 0; i < number_regions; i++) {
        model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Infected)}]  = 1;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Recovered)}] = 0;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Susceptible)}] =
            total_population_per_region -
            model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
                mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Infected)}] -
            model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
                mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Recovered)}];
    }

    model.parameters.set<mio::osirmobility::TimeInfected>(2);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.set<mio::osirmobility::ImpactCommuters>(1.);
    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 1.;
    model.parameters.get<mio::osirmobility::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.5});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(2), 0.8});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(2), mio::osirmobility::Region(0), 0.5});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(0), mio::osirmobility::Region(3), 1.0});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(3), 0.8});

    auto preprocess = run(filename, model);

    auto integrator = std::make_shared<mio::EulerIntegratorCore>();

    model.check_constraints();

    auto sir = simulate(t0, tmax, dt, model, integrator);

    bool print_to_terminal = true;

    sir.print_table();

    if (print_to_terminal) {

        std::vector<std::string> vars = {"S", "I", "R"};
        printf("Number of time points :%d\n", static_cast<int>(sir.get_num_time_points()));
        printf("People in\n");

        for (size_t k = 0; k < (size_t)mio::osirmobility::InfectionState::Count; k++) {
            double dummy = 0;

            for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
                printf("\t %s[%d]: %.2f", vars[k].c_str(), (int)i,
                       sir.get_last_value()[k + (size_t)mio::osirmobility::InfectionState::Count * (int)i]);
                dummy += sir.get_last_value()[k + (size_t)mio::osirmobility::InfectionState::Count * (int)i];
            }

            printf("\t %s_total: %.2f\n", vars[k].c_str(), dummy);
        }
    }
}