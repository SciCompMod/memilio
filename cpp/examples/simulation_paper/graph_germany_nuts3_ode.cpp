/*
* Copyright (C) 2020-2025 German Aerospace Center (DLR-SC)
*
* Authors: Carlotta Gerstein, Maximilian Betz
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
#include "memilio/io/cli.h"
#include "memilio/utils/base_dir.h"
#include "memilio/io/json_serializer.h"

#include "boost/filesystem.hpp"
#include <cstddef>
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Work = 0,
    Count,
};

enum class RunMode
{
    Load,
    Save,
};

enum class TestCase
{
    Open,
    KeepNPIFomInference,
    Lockdown,
    Dynamic,
};

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecir::Parameters<double>& params)
{

    //times
    // params.get<mio::osecir::TimeExposed<double>>()            = 3.3;
    // params.get<mio::osecir::TimeInfectedNoSymptoms<double>>() = 1.9;
    // params.get<mio::osecir::TimeInfectedSymptoms<double>>()   = 6.9;
    // params.get<mio::osecir::TimeInfectedSevere<double>>()     = 8;
    // params.get<mio::osecir::TimeInfectedCritical<double>>()   = 10.5;

    // //probabilities
    // params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()  = 0.08;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()    = 1;

    // params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()    = 0.207;
    // params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()         = 0.079;
    // params.get<mio::osecir::CriticalPerSevere<double>>()                 = 0.173;
    // params.get<mio::osecir::DeathsPerCritical<double>>()                 = 0.217;

    params.get<mio::osecir::TestAndTraceCapacity<double>>()                 = 0.0;


    params.set<mio::osecir::StartDay<double>>(0);
    params.set<mio::osecir::Seasonality<double>>(0.2);

    // params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()    = 0.0;
    // params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>() = 0.45;

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
    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(num_groups, num_groups, 7.95),
                           Eigen::MatrixXd::Zero(num_groups, num_groups));

    mio::unused(params);
    return mio::success();
}

// mio::IOResult<void> set_npis(mio::osecir::Parameters<double>& params, mio::Date start_date)
// {   
//     size_t num_groups = size_t(params.get_num_groups());
//     // npi from inference
//     auto start_damping = mio::Date(2020, 10, 8);
//     mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
//     contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.2), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)));


//     // change in npi at end of inference
//     start_damping = mio::Date(2020, 11, 30);
//     contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.0), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)));
//     //contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)));

//     params.get<mio::osecir::DynamicNPIsImplementationDelay<double>>() = 5;
//     auto dynamic_npi_dampings = std::vector<mio::DampingSampling<double>>();
//     dynamic_npi_dampings.push_back(mio::DampingSampling<double>(0.3, mio::DampingLevel(0),
//                                             mio::DampingType(0), mio::SimulationTime(mio::get_offset_in_days(start_damping, start_date)),
//                                             {0},  Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0)));

//     auto& dynamic_npis        = params.get<mio::osecir::DynamicNPIsInfectedSymptoms<double>>();
//     dynamic_npis.set_interval(mio::SimulationTime(1.0));
//     dynamic_npis.set_duration(mio::SimulationTime(14.0));
//     dynamic_npis.set_base_value(100'000);
//     dynamic_npis.set_threshold(200.0, dynamic_npi_dampings);

//     return mio::success();
// }

mio::IOResult<void> set_initial_compartments(mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>& params_graph, const fs::path& data_dir)
{
    // might be problematic if nodes in result are wrongly ordered
    BOOST_OUTCOME_TRY(auto&& simulation_results, mio::read_result((data_dir / "results_run0.h5").string()));

    for (size_t node_idx = 0; node_idx < params_graph.nodes().size(); ++node_idx) {
        auto& node = params_graph.nodes()[node_idx];
        // int id = node.id;
        auto pop_last_step = simulation_results[node_idx].get_groups().get_last_value();
        mio::unused(pop_last_step);
        for ( int infection_state = 0; infection_state != (int)mio::osecir::InfectionState::Count; infection_state++ ) {
            node.property.populations[{mio::AgeGroup(0), static_cast<mio::osecir::InfectionState>(infection_state)}] = double(pop_last_step[infection_state]);
            // node.property.populations[{mio::AgeGroup(0), static_cast<mio::osecir::InfectionState>(infection_state)}] = 3000;
        }
    }
    return mio::success();
}

mio::IOResult<void> set_sampled_parameters(mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>& params_graph, const fs::path& data_dir, TestCase test_case)
{

    BOOST_OUTCOME_TRY(auto&& parameter_list, mio::read_json((data_dir / "samples.json").string()));

    for (size_t node_idx = 0; node_idx < params_graph.nodes().size(); ++node_idx) {
        auto& node = params_graph.nodes()[node_idx];
        auto& params = node.property.parameters;

        // set "t_E", "t_C", "t_ISy", "t_ISev", "t_Cr", "mu_CR", "mu_IH", "mu_HU", "mu_UD", "transmission_prob"
        params.get<mio::osecir::TimeExposed<double>>()              = parameter_list["t_E"][0].asDouble();
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()   = 5.2 - parameter_list["t_E"][0].asDouble();
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()     = parameter_list["t_ISy"][0].asDouble();
        params.get<mio::osecir::TimeInfectedSevere<double>>()       = parameter_list["t_ISev"][0].asDouble();
        params.get<mio::osecir::TimeInfectedCritical<double>>()     = parameter_list["t_Cr"][0].asDouble();
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()   = parameter_list["mu_CR"][0].asDouble();
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()        = parameter_list["mu_IH"][0].asDouble();
        params.get<mio::osecir::CriticalPerSevere<double>>()                = parameter_list["mu_HU"][0].asDouble();
        params.get<mio::osecir::DeathsPerCritical<double>>()                = parameter_list["mu_UD"][0].asDouble();
        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>() = parameter_list["transmission_prob"][0].asDouble();

        int state = (int)mio::regions::get_state_id(node.id);
        // set "damping_values"
        // auto start_damping = mio::Date(2020, 10, 8);
        double damping_value;
        if (test_case == TestCase::Open){
            damping_value = 0;
        }
        else if (test_case == TestCase::KeepNPIFomInference){
            // t = 15, 30, 45
            damping_value = parameter_list["damping_values"][state][2].asDouble();
        }
        else if (test_case == TestCase::Lockdown) {
            damping_value = 1;
        }
        else if (test_case == TestCase::Dynamic) {
            damping_value = 0;
        }

        mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
        contact_matrix[0].add_damping(Eigen::MatrixXd::Constant(1, 1, damping_value), mio::SimulationTime<double>(0));
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
mio::IOResult<mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, const int num_days, const fs::path& data_dir, TestCase test_case)
{
    // global parameters
    const int num_groups = 1;
    auto end_date        = mio::offset_date_by_days(start_date, num_days);
    mio::osecir::Parameters<double> params(num_groups);
    params.get<mio::osecir::StartDay<double>>() = mio::get_day_in_year(start_date);

    // set default values independent from inference, may also be overwritten later
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(params));
    // BOOST_OUTCOME_TRY(set_npis(params, start_date));

    // create graph from params
    auto population_data_path =
        mio::path_join((data_dir / "Germany" / "pydata").string(), "county_current_population.json");

    BOOST_OUTCOME_TRY(auto&& node_ids, mio::get_node_ids(population_data_path, true, true));
    params.check_constraints();

    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> params_graph;
    for (size_t node_idx = 0; node_idx < node_ids.size(); ++node_idx) {
        params_graph.add_node(node_ids[node_idx], num_groups);
        params_graph.nodes()[node_idx].property.parameters = params;
    }

    // initialize model variables dependent on inference
    BOOST_OUTCOME_TRY(set_initial_compartments(params_graph, data_dir));
    BOOST_OUTCOME_TRY(set_sampled_parameters(params_graph, data_dir, test_case));
    
    // set edges of the graph
    const auto& read_function_edges = mio::read_mobility_plain;

    auto mobile_compartments     = {mio::osecir::InfectionState::Susceptible,
                                    mio::osecir::InfectionState::Exposed,
                                    mio::osecir::InfectionState::InfectedNoSymptoms,
                                    mio::osecir::InfectionState::InfectedSymptoms,
                                    mio::osecir::InfectionState::Recovered};

    const auto& set_edge_function =
        mio::set_edges<double, ContactLocation, mio::osecir::Model<double>, mio::MobilityParameters<double>,
                       mio::MobilityCoefficientGroup<double>, mio::osecir::InfectionState, decltype(read_function_edges)>;
    
    BOOST_OUTCOME_TRY(set_edge_function(mio::path_join(data_dir.string(), "Germany", "mobility", "commuter_mobility_2022.txt"), params_graph, mobile_compartments, 1,
                                        read_function_edges, std::vector<double>{1}, std::vector<std::vector<size_t>>{}));
    
    mio::unused(end_date, test_case);
    return mio::success(params_graph);
}

mio::IOResult<void> run(const int num_days_sim, mio::Date start_date, const std::string& data_dir,
                        const std::string& result_dir, RunMode mode, TestCase test_case)
{
    std::string test_case_name;
    if (test_case == TestCase::Open){
        test_case_name = "open";
    }
    else if (test_case == TestCase::KeepNPIFomInference){
        // t = 15, 30, 45
        test_case_name = "same";
    }
    else if (test_case == TestCase::Lockdown) {
        test_case_name = "lockdown";
    }
    else if (test_case == TestCase::Dynamic) {
        test_case_name = "dynamic";
    }

    //create or load graph
    auto save_dir = mio::path_join(result_dir, "graph_" + test_case_name);
    boost::filesystem::path sav_dir(save_dir);
    boost::filesystem::create_directories(sav_dir);

    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> params_graph;
    if (mode == RunMode::Save) {
        BOOST_OUTCOME_TRY(auto&& created, get_graph(start_date, num_days_sim, data_dir, test_case));
        BOOST_OUTCOME_TRY(mio::write_graph(created, save_dir));
        params_graph = created;
    }
    else {
        BOOST_OUTCOME_TRY(auto&& loaded, mio::read_graph<double, mio::osecir::Model<double>>(save_dir));
        params_graph = loaded;
    }

    // Copy graph for simulation
    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    mio::Graph<mio::SimulationNode<double, mio::osecir::Simulation<double>>, mio::MobilityEdge<double>> graph;
    for (auto& node : params_graph.nodes()) {
        graph.add_node(node.id, node.property);
    }
    for (auto& edge : params_graph.edges()) {
        graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    for (auto& node : graph.nodes()) {
        node.property.get_simulation().set_integrator_core(std::make_unique<mio::EulerIntegratorCore<double>>());
    }

    // auto sim = mio::make_mobility_sim<double>(0.0, 0.5, std::move(graph));
    auto sim = mio::make_no_mobility_sim<double>(0.0, std::move(graph));
    sim.advance(num_days_sim);

    std::vector<mio::TimeSeries<double>> results;
    results.reserve(sim.get_graph().nodes().size());
    std::transform(sim.get_graph().nodes().begin(), sim.get_graph().nodes().end(), std::back_inserter(results),
                   [](auto& n) {
                       return n.property.get_result();
                   });

    // std::vector<mio::TimeSeries<double>> results = mio::interpolate_simulation_result(sim.get_graph());
    // BOOST_OUTCOME_TRY(auto&& json_node, mio::serialize_json(sim.get_graph().nodes()[0].property.get_simulation().get_model()));
    // BOOST_OUTCOME_TRY(mio::write_json(mio::path_join(result_dir, "node0_" + test_case_name + ".json"), json_node));
    auto res = mio::save_result(results, county_ids, 1, mio::path_join(result_dir, "result_" + test_case_name + ".h5"));

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::err);
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                          .add<"DataDirectory">(mio::path_join(mio::base_dir(), "data/"))
                          .add<"ResultDirectory">(mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper/results"))
                          .add<"NumberSimulationDays">(60, {.alias = "n"})
                          .add<"RunMode">(RunMode::Save)
                          .add<"TestCase">(TestCase::Open)
                          .add<"StartDateDay">(30, {.alias="day"})
                          .add<"StartDateMonth">(11, {.alias="month"})
                          .add<"StartDateYear">(2020, {.alias="year"})
                          .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"DataDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }


    mio::Date start_date = mio::Date(cli_parameters.get<"StartDateYear">(), cli_parameters.get<"StartDateMonth">(), cli_parameters.get<"StartDateDay">());

    // std::string data_dir = "../../data";
    // std::string result_dir = "./results";

    boost::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    boost::filesystem::create_directories(res_dir);

    auto result =
        run(cli_parameters.get<"NumberSimulationDays">(), start_date, cli_parameters.get<"DataDirectory">(), cli_parameters.get<"ResultDirectory">(), 
            cli_parameters.get<"RunMode">(), cli_parameters.get<"TestCase">());
    return 0;
}