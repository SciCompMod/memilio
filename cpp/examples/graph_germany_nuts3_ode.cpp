/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors:
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
#include "ode_secir/parameters_io.h"
#include "ode_secir/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstddef>
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
// enum class ContactLocation
// {
//     Home = 0,
//     School,
//     Work,
//     Other,
//     Count,
// };
enum class ContactLocation
{
    Work,
    Count,
};

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecir::Parameters<double>& params)
{

    //times
    params.get<mio::osecir::TimeExposed<double>>()            = 3.3;
    params.get<mio::osecir::TimeInfectedNoSymptoms<double>>() = 1.9;
    params.get<mio::osecir::TimeInfectedSymptoms<double>>()   = 6.9;
    params.get<mio::osecir::TimeInfectedSevere<double>>()     = 8;
    params.get<mio::osecir::TimeInfectedCritical<double>>()   = 10.5;

    //probabilities
    params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()  = 0.08;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()    = 1;

    params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()    = 0.207;
    params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()         = 0.079;
    params.get<mio::osecir::CriticalPerSevere<double>>()                 = 0.173;
    params.get<mio::osecir::DeathsPerCritical<double>>()                 = 0.217;

    params.set<mio::osecir::StartDay>(0);
    params.set<mio::osecir::Seasonality<double>>(0.2);

    // params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()    = 0.0;
    // params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>() = 0.45;
    // params.get<mio::osecir::TestAndTraceCapacity<double>>()              = 35;

    return mio::success();
}

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(mio::osecir::Parameters<double>& params)
{
    size_t num_groups = size_t(params.get_num_groups());
    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, 7.95),
                           Eigen::MatrixXd::Zero(num_groups, num_groups));

    return mio::success();
}

mio::IOResult<void> set_npis(mio::osecir::Parameters<double>& params, mio::Date start_date)
{   
    size_t num_groups = size_t(params.get_num_groups());
    // npi from inference
    auto start_damping = mio::Date(2020, 10, 8);
    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.2), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)));


    // change in npi at end of inference
    start_damping = mio::Date(2020, 11, 30);
    contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.0), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)));
    //contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)));

    params.get<mio::osecir::DynamicNPIsImplementationDelay<double>>() = 5;
    auto dynamic_npi_dampings = std::vector<mio::DampingSampling<double>>();
    dynamic_npi_dampings.push_back(mio::DampingSampling<double>(0.3, mio::DampingLevel(0),
                                            mio::DampingType(0), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)),
                                            {0},  Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0)));

    auto& dynamic_npis        = params.get<mio::osecir::DynamicNPIsInfectedSymptoms<double>>();
    dynamic_npis.set_interval(mio::SimulationTime(1.0));
    dynamic_npis.set_duration(mio::SimulationTime(14.0));
    dynamic_npis.set_base_value(100'000);
    dynamic_npis.set_threshold(200.0, dynamic_npi_dampings);

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
mio::IOResult<mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, const int num_days, const fs::path& data_dir)
{
    // global parameters
    const int num_groups = 1;
    auto end_date        = mio::offset_date_by_days(start_date, num_days);
    mio::osecir::Parameters params(num_groups);
    params.get<mio::osecir::StartDay>() = mio::get_day_in_year(start_date);

    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(params));
    BOOST_OUTCOME_TRY(set_npis(params, start_date));

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> params_graph;
    const auto& read_function_nodes = mio::osecir::read_input_data_county<mio::osecir::Model<double>>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 2.5);
    auto scaling_factor_icu      = 1.0;
    auto mobile_compartments     = {mio::osecir::InfectionState::Susceptible,
                                    mio::osecir::InfectionState::Exposed,
                                    mio::osecir::InfectionState::InfectedNoSymptoms,
                                    mio::osecir::InfectionState::InfectedSymptoms};
    auto tnt_capacity_factor     = 0.0; // 7.5 / 100000.0;

    const auto& set_node_function =
        mio::set_nodes<mio::osecir::TestAndTraceCapacity<double>, mio::osecir::ContactPatterns<double>,
                       mio::osecir::Model<double>, mio::MobilityParameters<double>,
                       mio::osecir::Parameters<double>, decltype(read_function_nodes), decltype(node_id_function)>;
    const auto& set_edge_function =
        mio::set_edges<ContactLocation, mio::osecir::Model<double>, mio::MobilityParameters<double>,
                       mio::MobilityCoefficientGroup, mio::osecir::InfectionState, decltype(read_function_edges)>;

    auto population_data_path =
        mio::path_join((data_dir / "Germany" / "pydata").string(), "county_current_population.json");

    BOOST_OUTCOME_TRY(set_node_function(params, start_date, end_date, mio::path_join(data_dir.string(), "Germany", "pydata"), population_data_path, true,
                                        params_graph, read_function_nodes, node_id_function, scaling_factor_infected,
                                        scaling_factor_icu, tnt_capacity_factor, num_days, false, true));
    BOOST_OUTCOME_TRY(set_edge_function(mio::path_join(data_dir.string(), "Germany", "mobility", "commuter_mobility_2019.txt"), params_graph, mobile_compartments, 1,
                                        read_function_edges, std::vector<double>{1.0},
                                        std::vector<std::vector<size_t>>{}));
    return mio::success(params_graph);
}

enum class RunMode
{
    Load,
    Save,
};

mio::IOResult<void> run(RunMode mode, const int num_days_sim, mio::Date start_date, const std::string& data_dir,
                        const std::string& result_dir)
{
    //create or load graph
    auto save_dir = mio::path_join(result_dir, "graph");
    boost::filesystem::path sav_dir(save_dir);
    boost::filesystem::create_directories(sav_dir);

    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> params_graph;
    if (mode == RunMode::Save) {
        BOOST_OUTCOME_TRY(auto&& created, get_graph(start_date, num_days_sim, data_dir));
        BOOST_OUTCOME_TRY(mio::write_graph(created, save_dir));
        params_graph = created;
    }
    else {
        BOOST_OUTCOME_TRY(auto&& loaded, mio::read_graph<double, mio::osecir::Model<double>>(save_dir));
        params_graph = loaded;
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    mio::Graph<mio::SimulationNode<mio::osecir::Simulation<double>>, mio::MobilityEdge<double>> graph;

    for (auto& node : params_graph.nodes()) {
        BOOST_OUTCOME_TRY(set_npis(node.property.parameters, start_date));
        graph.add_node(node.id, node.property);
    }
    for (auto& edge : params_graph.edges()) {
        graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    auto sim = mio::make_mobility_sim<double>(0.0, 0.1, std::move(graph));
    sim.advance(num_days_sim);

    std::vector<mio::TimeSeries<double>> results = mio::interpolate_simulation_result(sim.get_graph());
    // BOOST_OUTCOME_TRY(mio::save_result(results, county_ids, 
    //     size_t(graph.nodes()[0].property.get_simulation().get_model().parameters.get_num_groups()), mio::path_join(result_dir, "result.h")));
    auto res = mio::save_result(results, county_ids, 1, "result5.h5");

    return mio::success();
}

int main(int argc, char** argv)
{

    mio::unused(argc, argv);
    mio::set_log_level(mio::LogLevel::warn);

    std::string data_dir = "../../data";

    mio::Date start_date             = mio::Date(2020, 10, 1);
    int num_days_sim                 = 120;

    std::string result_dir = "./results";

    boost::filesystem::path res_dir(result_dir);
    boost::filesystem::create_directories(res_dir);

    auto result =
        run(RunMode::Load, num_days_sim, start_date, data_dir, result_dir);
    return 0;
}