
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/model.h"
#include "ode_sir_mobility/parameters.h"
#include "ode_sir_mobility/regions.h"
#include "memilio/io/io.h"

mio::IOResult<std::vector<std::vector<std::vector<int>>>> read_path_mobility(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, mio::count_lines(filename));

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

template <typename FP = ScalarType>
mio::IOResult<void> preprocess(const std::string& filename, mio::osirmobility::Model<FP>& model)
{
    BOOST_OUTCOME_TRY(auto&& mobility_paths, read_path_mobility(filename));
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
                model.parameters.template get<mio::osirmobility::PathIntersections>()[{
                    mio::osirmobility::Region(i), mio::osirmobility::Region(j)}] = intersection_region;
            }
        }
    }
    return mio::success();
}

template <typename FP = ScalarType>
mio::IOResult<void> set_mobility_weights(const std::string& mobility_data, const std::string& trip_chains,
                                         mio::osirmobility::Model<FP>& model, size_t number_regions)
{
    BOOST_OUTCOME_TRY(preprocess(trip_chains, model));
    // mobility between nodes
    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter,
                      mio::read_mobility_plain(mobility_data + "mobility" + "commuter_migration_scaled.txt"));
    if (mobility_data_commuter.rows() != Eigen::Index(number_regions) ||
        mobility_data_commuter.cols() != Eigen::Index(number_regions)) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

    for (auto age = mio::AgeGroup(0); age < model.parameters.get_num_agegroups(); age++) {
        for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
            for (size_t county_idx_j = 0; county_idx_j < number_regions; ++county_idx_j) {
                //commuters
                auto population_i      = model.populations.get_group_total(mio::osirmobility::Region(county_idx_i));
                auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / population_i;
                if (commuter_coeff_ij > 4e-5) {
                    model.parameters.template get<mio::osirmobility::CommutingRatio>().push_back(
                        {mio::osirmobility::Region(county_idx_i), mio::osirmobility::Region(county_idx_j),
                         commuter_coeff_ij});
                }
            }
        }
    }
    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 50.;
    ScalarType dt   = 1;

    ScalarType number_regions              = 4;
    ScalarType number_age_groups           = 1;
    ScalarType total_population_per_region = 10;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& mobility_data   = "";
    const std::string& trip_chain_data = "";

    mio::osirmobility::Model<ScalarType> model(number_regions, number_age_groups);

    for (size_t i = 0; i < number_regions; i++) {
        model.populations[{mio::osirmobility::Region(i), mio::AgeGroup(0),
                           mio::osirmobility::InfectionState::Infected}]  = 1;
        model.populations[{mio::osirmobility::Region(i), mio::AgeGroup(0),
                           mio::osirmobility::InfectionState::Recovered}] = 0;
        model.populations[{mio::osirmobility::Region(i), mio::AgeGroup(0),
                           mio::osirmobility::InfectionState::Susceptible}] =
            total_population_per_region -
            model.populations[{mio::osirmobility::Region(i), mio::AgeGroup(0),
                               mio::osirmobility::InfectionState::Infected}] -
            model.populations[{mio::osirmobility::Region(i), mio::AgeGroup(0),
                               mio::osirmobility::InfectionState::Recovered}];
    }

    model.parameters.set<mio::osirmobility::TimeInfected<>>(2);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact<>>(0.04);
    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting<>>(1.);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::osirmobility::ContactPatterns<ScalarType>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(1.0);
    contact_matrix[0].add_damping(0.6, mio::SimulationTime(12.5));

    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.2});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(2), 0.6});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(2), mio::osirmobility::Region(0), 0.5});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(0), mio::osirmobility::Region(3), 1.0});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(3), 0.2});

    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(0),
                                                                  mio::osirmobility::Region(1)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(0),
                                                                  mio::osirmobility::Region(3)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(1),
                                                                  mio::osirmobility::Region(0)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(1),
                                                                  mio::osirmobility::Region(2)}] = {0};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(1),
                                                                  mio::osirmobility::Region(3)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(2),
                                                                  mio::osirmobility::Region(1)}] = {0};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(3),
                                                                  mio::osirmobility::Region(0)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(3),
                                                                  mio::osirmobility::Region(1)}] = {2};

    // auto result_preprocess = set_mobility_weights(mobility_data, trip_chain_data, model, number_regions);

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator =
        std::make_shared<mio::EulerIntegratorCore<ScalarType>>();

    model.check_constraints();

    auto sir = simulate(t0, tmax, dt, model, integrator);

    bool print_to_terminal = true;

    sir.print_table();

    if (print_to_terminal) {

        std::vector<std::string> vars = {"S", "I", "R"};
        printf("\n # t");
        for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
            for (size_t k = 0; k < (size_t)mio::osirmobility::InfectionState::Count; k++) {
                printf(" %s_%d", vars[k].c_str(), (int)i);
            }
        }

        auto num_points = static_cast<size_t>(sir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", sir.get_time(i));
            for (size_t k = 0; k < (size_t)model.parameters.get_num_regions(); k++) {
                for (size_t j = 0; j < (size_t)mio::osirmobility::InfectionState::Count; j++) {
                    printf(" %.14f", sir.get_value(i)[j + (size_t)mio::osirmobility::InfectionState::Count * (int)k]);
                }
            }
        }
        printf("\n");
    }
}
