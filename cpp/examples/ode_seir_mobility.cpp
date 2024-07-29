
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "models/ode_seir_mobility/infection_state.h"
#include "models/ode_seir_mobility/model.h"
#include "models/ode_seir_mobility/parameters.h"
#include "models/ode_seir_mobility/regions.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"

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
mio::IOResult<void> preprocess(const std::string& filename, mio::oseirmobility::Model<FP>& model)
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
            std::vector<mio::oseirmobility::Region> intersection_region(intersection_int.size(),
                                                                        mio::oseirmobility::Region(0));
            for (size_t k = 0; k < n_regions; k++) {
                if (k == i || k == j) {
                    continue;
                }
                std::sort(mobility_paths[k][j].begin(), mobility_paths[k][j].end());
                std::set_intersection(mobility_paths[i][j].begin(), mobility_paths[i][j].end(),
                                      mobility_paths[k][j].begin(), mobility_paths[k][j].end(),
                                      std::back_inserter(intersection_int));

                if (intersection_int.begin() != intersection_int.end()) {
                    intersection_region.push_back(mio::oseirmobility::Region(k));
                    intersection_int.pop_back();
                }
            }
            if (intersection_region.begin() != intersection_region.end()) {
                model.parameters.template get<mio::oseirmobility::PathIntersections>()[{
                    mio::oseirmobility::Region(i), mio::oseirmobility::Region(j)}] = intersection_region;
            }
        }
    }
    return mio::success();
}

template <typename FP = ScalarType>
mio::IOResult<void> set_mobility_weights(const std::string& mobility_data, const std::string& trip_chains,
                                         mio::oseirmobility::Model<FP>& model, size_t number_regions)
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
                auto population_i      = model.populations.get_group_total(mio::oseirmobility::Region(county_idx_i));
                auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / population_i;
                if (commuter_coeff_ij > 4e-5) {
                    model.parameters.template get<mio::oseirmobility::CommutingRatio>().push_back(
                        {mio::oseirmobility::Region(county_idx_i), mio::oseirmobility::Region(county_idx_j),
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
    ScalarType tmax = 10.;
    ScalarType dt   = 1;

    std::vector<int> region_ids  = {1001, 1002};
    ScalarType number_regions    = region_ids.size();
    ScalarType number_age_groups = 1;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& mobility_data   = "";
    const std::string& trip_chain_data = "";

    mio::oseirmobility::Model<ScalarType> model(number_regions, number_age_groups);
    model.populations[{mio::oseirmobility::Region(0), mio::AgeGroup(0), mio::oseirmobility::InfectionState::Exposed}] =
        10;
    model.populations[{mio::oseirmobility::Region(0), mio::AgeGroup(0),
                       mio::oseirmobility::InfectionState::Susceptible}] = 9990;
    model.populations[{mio::oseirmobility::Region(1), mio::AgeGroup(0), mio::oseirmobility::InfectionState::Exposed}] =
        0;
    model.populations[{mio::oseirmobility::Region(1), mio::AgeGroup(0),
                       mio::oseirmobility::InfectionState::Susceptible}] = 10000;

    model.parameters.set<mio::oseirmobility::TransmissionProbabilityOnContact<>>(1.);

    model.parameters.set<mio::oseirmobility::TimeExposed<>>(3.);
    model.parameters.set<mio::oseirmobility::TimeInfected<>>(5.);

    model.parameters.set<mio::oseirmobility::ImpactTransmissionDuringCommuting<>>(0.);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseirmobility::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.5, mio::SimulationTime(5));

    model.parameters.get<mio::oseirmobility::CommutingRatio>().push_back(
        {mio::oseirmobility::Region(0), mio::oseirmobility::Region(1), 0.01});
    model.parameters.get<mio::oseirmobility::CommutingRatio>().push_back(
        {mio::oseirmobility::Region(1), mio::oseirmobility::Region(0), 0.01});

    using DefaultIntegratorCore =
        mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>;

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<DefaultIntegratorCore>();

    model.check_constraints();

    auto result_from_sim = simulate(t0, tmax, dt, model, integrator);

    auto save_result_status =
        mio::save_result({result_from_sim}, region_ids, number_regions * number_age_groups, "ode_result.h5");

    // bool print_to_terminal = true;

    // result_from_sim.print_table();

    // if (print_to_terminal) {

    //     std::vector<std::string> vars = {"S", "E", "I", "R"};
    //     printf("\n # t");
    //     for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
    //         for (size_t k = 0; k < (size_t)mio::oseirmobility::InfectionState::Count; k++) {
    //             printf(" %s_%d", vars[k].c_str(), (int)i);
    //         }
    //     }

    //     auto num_points = static_cast<size_t>(sir.get_num_time_points());
    //     for (size_t i = 0; i < num_points; i++) {
    //         printf("\n%.14f ", sir.get_time(i));
    //         for (size_t k = 0; k < (size_t)model.parameters.get_num_regions(); k++) {
    //             for (size_t j = 0; j < (size_t)mio::oseirmobility::InfectionState::Count; j++) {
    //                 printf(" %.14f", sir.get_value(i)[j + (size_t)mio::oseirmobility::InfectionState::Count * (int)k]);
    //             }
    //         }
    //     }
    //     printf("\n");
    // }
}
