/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow, Martin Kühn
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "memilio/compartments/parameter_studies.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
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
void assign_uniform_distribution(mio::UncertainValue<double>& p, double min, double max)
{
    p = mio::UncertainValue<double>(0.5 * (max + min));
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
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
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
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       double min, double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<double>& params)
{
    //times
    // TimeExposed and TimeInfectedNoSymptoms are calculated as described in
    // Khailaie et al. (https://doi.org/10.1186/s12916-020-01884-4)
    // given SI_min = 3.935, SI_max = 4.6, INC = 5.2
    const double timeExposedMin            = 2.67;
    const double timeInfectedNoSymptomsMin = 1.2;

    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const double timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed<double>>(), timeExposedMin,
                                      timeExposedMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical<double>>(),
                                      timeInfectedCriticalMin, timeInfectedCriticalMin);

    //probabilities
    double fac_variant                                 = 1.4;
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};
    const double relativeTransmissionNoSymptomsMin     = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const double riskOfInfectionFromSymptomaticMin    = 0.0;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const double severePerInfectedSymptomsMin[]       = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const double criticalPerSevereMin[]               = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35};
    const double deathsPerCriticalMin[]               = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};

    const double reducExposedPartialImmunityMin                     = 0.75;
    const double reducExposedImprovedImmunityMin                    = 0.281;
    const double reducInfectedSymptomsPartialImmunityMin            = 0.6;
    const double reducInfectedSymptomsImprovedImmunityMin           = 0.193;
    const double reducInfectedSevereCriticalDeadPartialImmunityMin  = 0.05;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.041;

    double temp_reducTimeInfectedMild = 0.5;

    const double reducTimeInfectedMild = temp_reducTimeInfectedMild;

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMin);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMin);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMin);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMin);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild<double>>(),
                                      reducTimeInfectedMild, reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.1;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality<double>>(), seasonality_min, seasonality_min);

    // Delta specific parameter
    params.get<mio::osecirvvs::StartDayNewVariant>() = mio::get_day_in_year(mio::Date(2021, 6, 6));

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
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters<double>& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

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
mio::IOResult<mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, mio::Date end_date, const fs::path& data_dir)
{

    // global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters<double> params(num_age_groups);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    auto tnt_capacity_factor     = 1.43 / 100000.;
    auto mobile_compartments     = {mio::osecirvvs::InfectionState::SusceptibleNaive,
                                    mio::osecirvvs::InfectionState::ExposedNaive,
                                    mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive,
                                    mio::osecirvvs::InfectionState::InfectedSymptomsNaive,
                                    mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
                                    mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
                                    mio::osecirvvs::InfectionState::ExposedPartialImmunity,
                                    mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity,
                                    mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity,
                                    mio::osecirvvs::InfectionState::ExposedImprovedImmunity,
                                    mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity,
                                    mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity};

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    const auto& set_node_function =
        mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity<double>, mio::osecirvvs::ContactPatterns<double>,
                       mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::osecirvvs::Parameters<double>, decltype(read_function_nodes), decltype(node_id_function)>;
    const auto& set_edge_function =
        mio::set_edges<ContactLocation, mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::MobilityCoefficientGroup, mio::osecirvvs::InfectionState, decltype(read_function_edges)>;
    BOOST_OUTCOME_TRY(set_node_function(
        params, start_date, end_date, data_dir,
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json"), true,
        params_graph, read_function_nodes, node_id_function, scaling_factor_infected, scaling_factor_icu,
        tnt_capacity_factor, mio::get_offset_in_days(end_date, start_date), false, true));
    BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, mobile_compartments, contact_locations.size(),
                                        read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.},
                                        {}));

    return mio::success(params_graph);
}

/**
 * Different modes for running the parameter study.
 */
enum class RunMode
{
    Load,
    Save,
};

using ParamGraph      = mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>;
using SimulationGraph = mio::Graph<mio::SimulationNode<mio::osecirvvs::Simulation<>>, mio::MobilityEdge<double>>;

mio::IOResult<mio::GraphSimulation<SimulationGraph>>
create_graph_simulation(mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph,
                        double t0, double dt = 0.5)
{
    SimulationGraph sim_graph;

    for (auto&& node : params_graph.nodes()) {
        sim_graph.add_node(node.id, node.property, t0, dt);
    }
    for (auto&& edge : params_graph.edges()) {
        sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    return make_mobility_sim(t0, dt, std::move(sim_graph));
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
mio::IOResult<void> run(const fs::path& data_dir)
{

    const auto start_date   = mio::Date(2021, 6, 1);
    const auto num_days_sim = 5.0;
    const auto end_date     = mio::offset_date_by_days(start_date, int(std::ceil(num_days_sim)));

    //create graph
    BOOST_OUTCOME_TRY(auto&& params_graph, get_graph(start_date, end_date, data_dir));

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    // simulate the graph without using the parameter study
    // first transform the parameter graph into a simulation graph
    BOOST_OUTCOME_TRY(auto&& graph_simulation, create_graph_simulation(std::move(params_graph), 0.0));
    // if failed print error
    graph_simulation.advance(num_days_sim);

    // example how to get the results
    auto& results_graph = graph_simulation.get_graph();
    auto results_node_0 = results_graph.nodes()[0].property.get_result();
    std::cout << "last time point " << results_node_0.get_last_time() << "\n";

    return mio::success();
}

int main()
{

    mio::set_log_level(mio::LogLevel::warn);

    std::string data_dir = ""; // PATH TO DATA DIR

    auto result = run(data_dir);
    return 0;
}
