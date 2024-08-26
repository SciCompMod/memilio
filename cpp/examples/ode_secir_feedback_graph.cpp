#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/config.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secir/parameters_io.h"
#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

/**
 * different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
};

/**
 * different level of NPI, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
};

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue& p, double min, double max)
{
    p = mio::UncertainValue(0.5 * (max + min));
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min[i] and max[i] as a value and UNIFORM(min[i], max[i]) as a distribution for
 * each element i of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution for each element of array.
 * @param max minimum of distribution for each element of array.
 */
template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution to every element of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double min,
                                       double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecir::Parameters& params)
{
    //times
    // TimeExposed and TimeInfectedNoSymptoms are calculated as described in
    // Khailaie et al. (https://doi.org/10.1186/s12916-020-01884-4)
    // given SI_min = 3.935, SI_max = 4.6, INC = 5.2
    const double timeExposedMin = 2.67;
    // const double timeExposedMax            = 4.;
    const double timeInfectedNoSymptomsMin = 1.2;
    // const double timeInfectedNoSymptomsMax = 2.53;

    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    // const double timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const double reduc_time_severe = 0.5;
    double timeInfectedSevereMin[] = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    // double timeInfectedSevereMax[]         = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    // multiply timeInfectedSevereMin and timeInfectedSevereMax by reduc_time_severe
    for (size_t i = 0; i < 6; i++) {
        timeInfectedSevereMin[i] *= reduc_time_severe;
        // timeInfectedSevereMax[i] *= reduc_time_severe;
    }
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    // const double timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

    array_assign_uniform_distribution(params.get<mio::osecir::TimeExposed>(), timeExposedMin, timeExposedMin);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedNoSymptoms>(), timeInfectedNoSymptomsMin,
                                      timeInfectedNoSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedSymptoms>(), timeInfectedSymptomsMin,
                                      timeInfectedSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedSevere>(), timeInfectedSevereMin,
                                      timeInfectedSevereMin);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedCritical>(), timeInfectedCriticalMin,
                                      timeInfectedCriticalMin);

    //probabilities
    double fact = 0.8;
    // const double transmissionProbabilityOnContactMin[] = {fact * 0.02, fact * 0.05, fact * 0.05,
    //                                                       fact * 0.05, fact * 0.08, fact * 0.15};
    const double transmissionProbabilityOnContactMax[] = {fact * 0.04, fact * 0.07, fact * 0.07,
                                                          fact * 0.07, fact * 0.10, fact * 0.20};
    const double relativeTransmissionNoSymptomsMin     = 1;
    // const double relativeTransmissionNoSymptomsMax     = 1;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    // const double riskOfInfectionFromSymptomaticMin = 0.1;
    const double riskOfInfectionFromSymptomaticMax = 0.3;
    // const double maxRiskOfInfectionFromSymptomaticMin = 0.3;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;
    // const double recoveredPerInfectedNoSymptomsMin[] = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const double recoveredPerInfectedNoSymptomsMax[] = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    // const double severePerInfectedSymptomsMin[] = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const double severePerInfectedSymptomsMax[] = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double fact_severe                    = 1.5;
    // const double criticalPerSevereMin[] = {fact_severe * 0.05, fact_severe * 0.05, fact_severe * 0.05,
    //                                        fact_severe * 0.10, fact_severe * 0.25, fact_severe * 0.35};
    const double criticalPerSevereMax[] = {fact_severe * 0.10, fact_severe * 0.10, fact_severe * 0.10,
                                           fact_severe * 0.20, fact_severe * 0.35, fact_severe * 0.45};
    const double deathsPerCriticalMin[] = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};
    // const double deathsPerCriticalMax[]               = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    array_assign_uniform_distribution(params.get<mio::osecir::TransmissionProbabilityOnContact>(),
                                      transmissionProbabilityOnContactMax, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RelativeTransmissionNoSymptoms>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecir::RiskOfInfectionFromSymptomatic>(),
                                      riskOfInfectionFromSymptomaticMax, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>(),
                                      maxRiskOfInfectionFromSymptomaticMax, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>(),
                                      recoveredPerInfectedNoSymptomsMax, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::SeverePerInfectedSymptoms>(),
                                      severePerInfectedSymptomsMax, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::CriticalPerSevere>(), criticalPerSevereMax,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecir::DeathsPerCritical>(), deathsPerCriticalMin,
                                      deathsPerCriticalMin);

    //sasonality
    // const double seasonality_min = 0.1;
    // const double seasonality_max = 0.3;

    // assign_uniform_distribution(params.get<mio::osecir::Seasonality>(), seasonality_min, seasonality_max);

    return mio::success();
}

/**
 * Set feedback parameters
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_feedback_parameters(mio::osecir::Parameters& params, ScalarType min, ScalarType max,
                                            ScalarType icu_capacity, ScalarType fact_regional)
{
    params.get<mio::osecir::ICUCapacity>()            = icu_capacity;
    params.get<mio::osecir::CutOffGamma>()            = 45;
    params.get<mio::osecir::EpsilonContacts>()        = 0.1;
    params.get<mio::osecir::BlendingFactorLocal>()    = 0.3;
    params.get<mio::osecir::BlendingFactorRegional>() = fact_regional;
    params.get<mio::osecir::ContactReductionMin>()    = {min, min, min, min};
    params.get<mio::osecir::ContactReductionMax>()    = {max, max, max, max};
    auto& icu_occupancy                               = params.get<mio::osecir::ICUOccupancyLocal>();
    Eigen::VectorXd icu_day                           = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 0.0);
    for (int t = -50; t <= 0; ++t) {
        icu_occupancy.add_time_point(t, icu_day);
    }

    return mio::success();
}

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecir::Parameters& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        // BOOST_OUTCOME_TRY(auto&& minimum,
        //                   mio::read_mobility_plain(
        //                       (data_dir / "contacts" / ("minimum_" + contact_location.second + ".txt")).string()));
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrices);

    return mio::success();
}

void set_state_ids(mio::Graph<mio::osecir::Model, mio::MigrationParameters>& graph)
{
    for (auto& node : graph.nodes()) {
        int id                    = node.id;
        std::string county_id_str = std::to_string(id);
        if (county_id_str.length() == 4) {
            // four digit numbers: first digit is the State-ID
            id = county_id_str[0] - '0';
        }
        else if (county_id_str.length() == 5) {
            // five digit numbers: first two digits are the State-ID
            id = std::stoi(county_id_str.substr(0, 2));
        }
        else {
            throw std::invalid_argument("Invalid county_id length. Must be 4 or 5 digits long.");
        }
        node.property.parameters.get<mio::osecir::StateID>() = id;
    }
}

/**
 * Set NPIs.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param params Object that the NPIs will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_npis(mio::osecir::Parameters& params, const std::string& mode,
                             const ScalarType reduc_fac_location)
{
    if (std::strcmp(mode.c_str(), "ClassicDamping") == 0) {
        auto& contacts         = params.get<mio::osecir::ContactPatterns>();
        auto& contact_dampings = contacts.get_dampings();

        const size_t locations       = static_cast<size_t>(ContactLocation::Count);
        const auto group_weights_all = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0);

        auto contact_reduction = [=](auto val, size_t location) {
            auto v = mio::UncertainValue(val);
            return mio::DampingSampling(v, mio::DampingLevel(int(0)), mio::DampingType(int(0)), mio::SimulationTime(0),
                                        {location}, group_weights_all);
        };

        for (size_t loc = 0; loc < locations; ++loc) {
            contact_dampings.push_back(contact_reduction(reduc_fac_location, loc));
        }
        contacts.make_matrix();
    }

    return mio::success();
}

/**
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<mio::Graph<mio::osecir::Model, mio::MigrationParameters>>
get_graph(mio::Date start_date, mio::Date end_date, const fs::path& data_dir, const std::string& mode, ScalarType kmin,
          ScalarType kmax, ScalarType icu_capacity, ScalarType fact_regional)
{
    const auto start_day = mio::get_day_in_year(start_date);

    // global parameters
    const int num_age_groups = 6;
    mio::osecir::Parameters params(num_age_groups);
    params.get<mio::osecir::StartDay>() = start_day;
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    if (std::strcmp(mode.c_str(), "FeedbackDamping") == 0) {
        BOOST_OUTCOME_TRY(set_feedback_parameters(params, kmin, kmax, icu_capacity, fact_regional));
    }
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_npis(params, mode, kmin));

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1);
    auto scaling_factor_icu      = 1.0;
    auto tnt_capacity_factor     = 7.5 / 100000.;
    auto migrating_compartments  = {mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed,
                                   mio::osecir::InfectionState::InfectedNoSymptoms,
                                   mio::osecir::InfectionState::InfectedSymptoms,
                                   mio::osecir::InfectionState::Recovered};

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecir::Model, mio::MigrationParameters> params_graph;
    const auto& read_function_nodes = mio::osecir::read_input_data_county<mio::osecir::Model>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    const auto& set_node_function =
        mio::set_nodes<mio::osecir::TestAndTraceCapacity, mio::osecir::ContactPatterns, mio::osecir::Model,
                       mio::MigrationParameters, mio::osecir::Parameters, decltype(read_function_nodes),
                       decltype(node_id_function)>;
    const auto& set_edge_function =
        mio::set_edges<ContactLocation, mio::osecir::Model, mio::MigrationParameters, mio::MigrationCoefficientGroup,
                       mio::osecir::InfectionState, decltype(read_function_edges)>;
    BOOST_OUTCOME_TRY(
        set_node_function(params, start_date, end_date, data_dir,
                          mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json"),
                          true, params_graph, read_function_nodes, node_id_function, scaling_factor_infected,
                          scaling_factor_icu, tnt_capacity_factor, 0, false, true));
    BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, migrating_compartments, contact_locations.size(),
                                        read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}));

    set_state_ids(params_graph);
    return mio::success(params_graph);
}

mio::TimeSeries<ScalarType>
get_risk(mio::osecir::FeedbackSimulation<mio::FlowSimulation<mio::osecir::Model>>& simulation)
{
    auto& risk = simulation.get_perceived_risk();
    return mio::interpolate_simulation_result(risk);
}

mio::TimeSeries<ScalarType> get_risk(mio::osecir::Simulation<mio::FlowSimulation<mio::osecir::Model>>& /*simulation*/)
{
    return mio::TimeSeries<ScalarType>(6);
}

template <typename ParameterStudy>
mio::IOResult<void> run_parameter_study(ParameterStudy parameter_study, std::vector<int> county_ids,
                                        const fs::path& result_dir, const ScalarType num_days_sim,
                                        const std::string mode)
{

    // parameter_study.get_rng().seed(
    //    {114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203}); //set seeds, e.g., for debugging
    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : parameter_study.get_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }

    auto ensemble = parameter_study.run(
        [](auto&& graph) {
            return draw_sample(graph);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);

            auto params = std::vector<mio::osecir::Model>{};
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                           [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            auto flows = std::vector<mio::TimeSeries<ScalarType>>{};
            flows.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(flows),
                           [](auto&& node) {
                               auto& flow_node         = node.property.get_simulation().get_flows();
                               auto interpolated_flows = mio::interpolate_simulation_result(flow_node);
                               return interpolated_flows;
                           });

            // save the perceived risk for each run
            auto risks = std::vector<mio::TimeSeries<ScalarType>>{};
            if (std::strcmp(mode.c_str(), "FeedbackDamping") == 0) {
                risks.reserve(results_graph.nodes().size());
                std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(risks),
                               [](auto&& node) {
                                   auto interpolated_risk = get_risk(node.property.get_simulation());
                                   return mio::interpolate_simulation_result(interpolated_risk);
                               });
            }

            const auto num_groups =
                (size_t)results_graph.nodes()[0].property.get_simulation().get_model().parameters.get_num_groups();
            auto r0 = std::vector<mio::TimeSeries<ScalarType>>{};
            r0.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(r0),
                           [num_groups](auto&& node) {
                               auto r0_node   = mio::TimeSeries<ScalarType>(num_groups);
                               auto r0_values = mio::osecir::get_reproduction_numbers(node.property.get_simulation());
                               for (int i = 0; i < node.property.get_simulation().get_result().get_num_time_points();
                                    i++) {
                                   r0_node.add_time_point(node.property.get_simulation().get_result().get_time(i),
                                                          Eigen::VectorXd::Constant(num_groups, r0_values[i]));
                               }
                               auto interpolated_r0 = mio::interpolate_simulation_result(r0_node);
                               return mio::interpolate_simulation_result(interpolated_r0);
                           });

            auto contacts = std::vector<mio::TimeSeries<ScalarType>>{};
            contacts.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(contacts),
                           [num_groups, num_days_sim](auto&& node) {
                               Eigen::VectorXd contact_rate_ages(num_groups);
                               auto contacts_node = mio::TimeSeries<ScalarType>(num_groups);
                               mio::ContactMatrixGroup const& contact_matrix =
                                   node.property.get_simulation()
                                       .get_model()
                                       .parameters.template get<mio::osecir::ContactPatterns>();
                               for (auto t_sim = 0.0; t_sim < num_days_sim; t_sim++) {
                                   auto t_sim_day = mio::SimulationTime(t_sim);
                                   for (size_t i = 0; i < num_groups; ++i) {
                                       for (size_t j = 0; j < num_groups; ++j) {
                                           contact_rate_ages(i) += contact_matrix.get_matrix_at(t_sim_day)(
                                               static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                                       }
                                   }
                                   contacts_node.add_time_point(t_sim, contact_rate_ages);
                               }
                               auto interpolated_contacts = mio::interpolate_simulation_result(contacts_node);
                               return interpolated_contacts;
                           });

            // for (size_t i = 0; i < results_graph.nodes().size(); i++) {
            //     auto& node = results_graph.nodes()[i];
            //     std::stringstream output1;
            //     node.property.get_simulation()
            //         .get_model_ptr()
            //         ->parameters.template get<mio::osecir::ICUOccupancyLocal>()
            //         .print_table({}, 12, 15, output1);

            //     std::string results_dir = "/localdata1/code_2024/memilio/test";

            //     std::ofstream file1(results_dir + "/icu_occupancy_local_" + std::to_string(county_ids[i]) + ".txt");
            //     file1 << output1.str();
            //     file1.close();
            // }

            std::cout << "run " << run_idx << " done" << std::endl;

            return std::make_tuple(std::move(interpolated_result), std::move(params), std::move(flows),
                                   std::move(risks), std::move(r0), std::move(contacts));
        });
    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecir::Model>>{};
        ensemble_params.reserve(ensemble.size());
        auto ensemble_flows = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_flows.reserve(ensemble.size());
        auto ensemble_risks = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_risks.reserve(ensemble.size());
        auto ensemble_r0 = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_r0.reserve(ensemble.size());
        auto ensemble_contacts = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_contacts.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(std::get<0>(run)));
            ensemble_params.emplace_back(std::move(std::get<1>(run)));
            ensemble_flows.emplace_back(std::move(std::get<2>(run)));
            ensemble_risks.emplace_back(std::move(std::get<3>(run)));
            ensemble_r0.emplace_back(std::move(std::get<4>(run)));
            ensemble_contacts.emplace_back(std::move(std::get<5>(run)));
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, county_ids, result_dir, false));

        auto result_dir_run_flows = result_dir / "flows";
        if (mio::mpi::is_root()) {
            boost::filesystem::create_directories(result_dir_run_flows);
            printf("Saving Flow results to \"%s\".\n", result_dir_run_flows.c_str());
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_flows, ensemble_params, county_ids, result_dir_run_flows, false));

        if (std::strcmp(mode.c_str(), "FeedbackDamping") == 0) {
            auto result_dir_risk = result_dir / "risk";
            if (mio::mpi::is_root()) {
                boost::filesystem::create_directories(result_dir_risk);
                printf("Saving Risk results to \"%s\".\n", result_dir_risk.c_str());
            }
            BOOST_OUTCOME_TRY(save_results(ensemble_risks, ensemble_params, county_ids, result_dir_risk, false));
        }

        auto result_dir_r0 = result_dir / "r0";
        if (mio::mpi::is_root()) {
            boost::filesystem::create_directories(result_dir_r0);
            printf("Saving R0 results to \"%s\".\n", result_dir_r0.c_str());
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_r0, ensemble_params, county_ids, result_dir_r0, false));

        auto result_dir_contacts = result_dir / "contacts";
        if (mio::mpi::is_root()) {
            boost::filesystem::create_directories(result_dir_contacts);
            printf("Saving Contact results to \"%s\".\n", result_dir_contacts.c_str());
        }
        BOOST_OUTCOME_TRY(save_results(ensemble_contacts, ensemble_params, county_ids, result_dir_contacts, false));
    }

    return mio::success();
}

/**
 * Run the parameter study.
 * Load a previously stored graph or create a new one from data. 
 * The graph is the input for the parameter study.
 * A newly created graph is saved and can be reused.
 * @param mode Mode for running the parameter study.
 * @param data_dir data directory. Not used if mode is RunMode::Load.
 * @param save_dir directory where the graph is loaded from if mode is RunMOde::Load or save to if mode is RunMode::Save.
 * @param result_dir directory where all results of the parameter study will be stored.
 * @param save_single_runs [Default: true] Defines if single run results are written to the disk.
 * @returns any io error that occurs during reading or writing of files.
 */
mio::IOResult<void> run(const fs::path& data_dir, const fs::path& result_dir)
{
    const auto start_date = mio::Date(2020, 10, 1);

    const auto num_days_sim = 200.0;
    const auto end_date     = mio::offset_date_by_days(start_date, int(std::ceil(num_days_sim)));
    const auto num_runs     = 1;

    auto const modes = {"FeedbackDamping", "ClassicDamping"}; //, "FeedbackDamping"};

    auto min_values = std::vector<ScalarType>{0.0};

    // auto max_values = std::vector<ScalarType>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    auto max_values = std::vector<ScalarType>{0.02, 0.04, 0.06, 0.08, 0.12, 0.14, 0.16, 0.18, 0.22, 0.24,
                                              0.26, 0.28, 0.32, 0.34, 0.36, 0.38, 0.42, 0.44, 0.46, 0.48};

    auto icu_capacities = std::vector<ScalarType>{9.0};

    // auto fact_regional = std::vector<ScalarType>{0.6, 0.7};
    auto fact_regional = std::vector<ScalarType>{0.02, 0.04, 0.06, 0.08, 0.12, 0.14, 0.16, 0.18,
                                                 0.22, 0.24, 0.26, 0.28, 0.32, 0.34, 0.36, 0.38};
    // const size_t county_id_infected = 3241;

    //create or load graph
    for (auto icu_cap : icu_capacities) {
        for (auto fact_r : fact_regional) {
            for (auto mode : modes) {
                for (size_t min_indx = 0; min_indx < min_values.size(); min_indx++) {
                    for (size_t max_indx = min_indx; max_indx < max_values.size(); max_indx++) {
                        if (std::strcmp(mode, "ClassicDamping") == 0 && min_values[min_indx] != max_values[max_indx]) {
                            continue;
                        }
                        auto& kmin = min_values[min_indx];
                        auto& kmax = max_values[max_indx];

                        mio::Graph<mio::osecir::Model, mio::MigrationParameters> params_graph;
                        auto result_dir_mode = result_dir / ("ICUCap_" + std::to_string(icu_cap)) /
                                               ("BlendingFactorRegional_" + std::to_string(fact_r)) /
                                               ("kmin_" + std::to_string(kmin) + "_kmax_" + std::to_string(kmax)) /
                                               boost::filesystem::path(mode);
                        if (std::strcmp(mode, "ClassicDamping") == 0) {
                            result_dir_mode = result_dir / ("ICUCap_" + std::to_string(icu_cap)) /
                                              ("BlendingFactorRegional_" + std::to_string(fact_r)) /
                                              ("fixed_damping_kmin_" + std::to_string(kmin)) / "ClassicDamping";
                        }
                        // create directory for results
                        if (mio::mpi::is_root()) {
                            boost::filesystem::create_directories(result_dir_mode);
                        }

                        BOOST_OUTCOME_TRY(auto&& created,
                                          get_graph(start_date, end_date, data_dir, mode, kmin, kmax, icu_cap, fact_r));
                        params_graph = created;

                        std::vector<int> county_ids(params_graph.nodes().size());
                        std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(),
                                       [](auto& n) {
                                           return n.id;
                                       });

                        //run parameter study
                        if (std::strcmp(mode, "ClassicDamping") == 0) {
                            auto parameter_study_sim =
                                mio::ParameterStudy<mio::osecir::Simulation<mio::FlowSimulation<mio::osecir::Model>>>{
                                    params_graph, 0.0, num_days_sim, 0.5, size_t(num_runs)};
                            BOOST_OUTCOME_TRY(run_parameter_study(parameter_study_sim, county_ids, result_dir_mode,
                                                                  num_days_sim, mode));
                        }
                        else {
                            auto parameter_study_feedback_sim = mio::ParameterStudy<
                                mio::osecir::FeedbackSimulation<mio::FlowSimulation<mio::osecir::Model>>>{
                                params_graph, 0.0, num_days_sim, 0.5, size_t(num_runs)};
                            BOOST_OUTCOME_TRY(run_parameter_study(parameter_study_feedback_sim, county_ids,
                                                                  result_dir_mode, num_days_sim, mode));
                        }
                    }
                }
            }
        }
    }
    return mio::success();
}
int main()
{
    mio::set_log_level(mio::LogLevel::err);
    mio::mpi::init();

    std::string data_dir   = "/localdata1/code_2024/memilio/data";
    std::string result_dir = "/localdata1/code_2024/memilio/results";

    auto result = run(data_dir, result_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}
