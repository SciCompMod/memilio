/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker, Daniel Abele, Martin Kühn
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

#include "memilio/utils/logging.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/mobility/graph.h"
#include "ode_secir/parameters_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/data/analyze_result.h"

namespace fs = boost::filesystem;

/**
 * Different modes for running the parameter study.
 */
enum class RunMode
{
    Load,
    Save,
};

/**
 * Set epidemiological parameters of Sars-CoV-2.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecir::Parameters<double>& params)
{
    params.get<mio::osecir::TimeExposed<double>>()            = 3.2;
    params.get<mio::osecir::TimeInfectedNoSymptoms<double>>() = 2.0;
    params.get<mio::osecir::TimeInfectedSymptoms<double>>()   = 5.8;
    params.get<mio::osecir::TimeInfectedSevere<double>>()     = 9.5;
    params.get<mio::osecir::TimeInfectedCritical<double>>()   = 7.1;

    params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()  = 0.05;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()    = 0.7;
    params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()    = 0.09;
    params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()    = 0.25;
    params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>() = 0.45;
    params.get<mio::osecir::TestAndTraceCapacity<double>>()              = 35;
    params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()         = 0.2;
    params.get<mio::osecir::CriticalPerSevere<double>>()                 = 0.25;
    params.get<mio::osecir::DeathsPerCritical<double>>()                 = 0.3;

    params.set<mio::osecir::Seasonality<double>>(0.2);

    return mio::success();
}

mio::IOResult<void> set_nodes(mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>& params_graph,
                              const mio::osecir::Parameters<double>& params, const fs::path& data_dir,
                              mio::Date start_date)
{
    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;

    //As we have only one age group in the json file the age group name is set to 'All'
    mio::ConfirmedCasesDataEntry::age_group_names = {"All"};

    mio::PopulationDataEntry::age_group_names = {
        "0-4 years",   "5-9 years",   "10-14 years", "15-19 years", "20-24 years", "25-29 years", "30-34 years",
        "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years", "65-69 years",
        "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90+ years"};

    //read node ids
    bool is_node_for_county         = false;
    bool interpolate_rki_age_groups = false;
    BOOST_OUTCOME_TRY(auto&& node_ids, mio::get_node_ids(mio::path_join((data_dir).string(), "population_data.json"),
                                                         is_node_for_county, interpolate_rki_age_groups));
    std::vector<mio::osecir::Model<double>> nodes(node_ids.size(),
                                                  mio::osecir::Model<double>(int(size_t(params.get_num_groups()))));
    //set parameters for every node
    for (auto& node : nodes) {
        node.parameters = params;
    }
    int num_days = 90;
    auto read = mio::osecir::read_input_data(nodes, start_date, node_ids, scaling_factor_infected, scaling_factor_icu,
                                             data_dir.string(), num_days, true);
    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        params_graph.add_node(node_ids[node_idx], nodes[node_idx]);
    }
    return mio::success();
}

mio::IOResult<void> set_edges(mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>& params_graph,
                              const fs::path& mobility_data_dir)
{
    auto mobile_compartments = {mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed,
                                mio::osecir::InfectionState::InfectedNoSymptoms,
                                mio::osecir::InfectionState::InfectedSymptoms, mio::osecir::InfectionState::Recovered};
    //mobility matrix has to be provided by the user as input and should have shape num_nodes x num_nodes
    BOOST_OUTCOME_TRY(auto&& mobility_data,
                      mio::read_mobility_plain(mio::path_join((mobility_data_dir).string(), "mobility_matrix.txt")));
    if (mobility_data.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data.cols() != Eigen::Index(params_graph.nodes().size())) {
        return mio::failure(mio::StatusCode::InvalidValue, "Mobility matrices do not have the correct size.");
    }

    for (size_t node_i = 0; node_i < params_graph.nodes().size(); ++node_i) {
        for (size_t node_j = 0; node_j < params_graph.nodes().size(); ++node_j) {
            auto& populations    = params_graph.nodes()[node_i].property.populations;
            auto mobility_coeffs = mio::MobilityCoefficientGroup<double>(1, populations.numel());

            auto coeff =
                mobility_data(node_i, node_j) < 1 ? 0 : mobility_data(node_i, node_j) / populations.get_total();

            for (auto age = mio::AgeGroup(0); age < populations.template size<mio::AgeGroup>(); ++age) {
                for (auto compartment : mobile_compartments) {
                    auto coeff_idx                               = populations.get_flat_index({age, compartment});
                    mobility_coeffs[0].get_baseline()[coeff_idx] = coeff;
                }
                if (coeff > 1e-5) {
                    params_graph.add_edge(node_i, node_j, std::move(mobility_coeffs));
                }
            }
        }
    }
    return mio::success();
}

/**
 *
*/
mio::IOResult<mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, const fs::path& data_dir)
{
    const int num_age_groups = 1;
    mio::osecir::Parameters<double> params(num_age_groups);
    params.get<mio::osecir::StartDay<double>>() = mio::get_day_in_year(start_date);

    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    double contact_freq                             = 10;
    mio::ContactMatrixGroup<double>& contact_matrix = params.template get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, contact_freq));
    contact_matrix.add_damping(0.7, mio::SimulationTime<double>(30));
    contact_matrix.add_damping(0.1, mio::SimulationTime<double>(50));

    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> params_graph;
    BOOST_OUTCOME_TRY(set_nodes(params_graph, params, data_dir, start_date));
    BOOST_OUTCOME_TRY(set_edges(params_graph, data_dir));

    return mio::success(params_graph);
}

/**
 * Run a graph simulation for Munich.
 * Load a previously stored graph or create a new one from data.
 * A newly created graph is saved and can be reused.
 * @param mode Mode for running the parameter study.
 * @param data_dir data directory. Not used if mode is RunMode::Load.
 * @param save_dir directory where the graph is loaded from if mode is RunMOde::Load or save to if mode is RunMode::Save.
 * @param result_dir directory where all results of the parameter study will be stored.
 * @returns any io error that occurs during reading or writing of files.
*/
mio::IOResult<void> run(RunMode mode, const fs::path& data_dir, const fs::path& save_dir, const fs::path& result_dir)
{
    const mio::Date start_date = mio::Date(2022, 1, 1);
    const auto num_days        = 90.0;

    //create or load graph
    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> params_graph;

    if (mode == RunMode::Save) {
        BOOST_OUTCOME_TRY(auto&& created_graph, get_graph(start_date, data_dir));
        BOOST_OUTCOME_TRY(mio::write_graph(created_graph, save_dir.string()));
        params_graph = created_graph;
    }
    else {
        BOOST_OUTCOME_TRY(auto&& loaded_graph, mio::read_graph<double, mio::osecir::Model<double>>(save_dir.string()));
        params_graph = loaded_graph;
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    //create simulation graph
    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecir::Model<double>>>,
               mio::MobilityEdge<double>>
        sim_graph;

    for (auto&& node : params_graph.nodes()) {
        sim_graph.add_node(node.id, node.property, 0.0);
    }
    for (auto&& edge : params_graph.edges()) {
        sim_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge.property);
    }

    auto sim = mio::make_mobility_sim<double>(0.0, 0.5, std::move(sim_graph));
    sim.advance(num_days);

    auto params = std::vector<mio::osecir::Model<double>>{};
    params.reserve(sim.get_graph().nodes().size());
    std::transform(sim.get_graph().nodes().begin(), sim.get_graph().nodes().end(), std::back_inserter(params),
                   [](auto&& node) {
                       return node.property.get_simulation().get_model();
                   });

    auto interpolated_result = mio::interpolate_simulation_result(std::move(sim).get_graph());

    auto save_single_run_result = save_result_with_params(interpolated_result, params, county_ids, result_dir, 0);

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);

    std::cout << "Current path is " << fs::current_path() << '\n';

    RunMode mode;
    std::string save_dir;
    std::string data_dir;
    std::string result_dir;
    bool save_single_runs = true;
    if (argc == 1) {
        mode       = RunMode::Save;
        save_dir   = "../../../munich_simulation/results";
        result_dir = "../../../munich_simulation/results";
        data_dir   = "../../../munich_simulation/data";
    }
    else if (argc == 4) {
        mode       = RunMode::Load;
        save_dir   = argv[1];
        result_dir = argv[2];
        data_dir   = "";
        printf("Loading graph from \"%s\".\n", save_dir.c_str());
        printf("Exporting single run results and parameters: %d.\n", (int)save_single_runs);
    }
    else if (argc == 5) {
        mode       = RunMode::Save;
        data_dir   = argv[1];
        save_dir   = argv[2];
        result_dir = argv[3];
        if (atoi(argv[4]) == 0) {
            save_single_runs = false;
        }
        printf("Reading data from \"%s\", saving graph to \"%s\".\n", data_dir.c_str(), save_dir.c_str());
        printf("Exporting single run results and parameters: %d.\n", (int)save_single_runs);
    }
    else {
        printf("Usage:\n");
        printf("munich_graph_sim <data_dir> <save_dir> <result_dir>");
        printf("\tMake graph with data from <data_dir> and save at <save_dir>, then run the simulation.\n");
        printf("\tStore the results in <result_dir>\n");
        printf("2021_vaccination_delta <load_dir> <result_dir>\n");
        printf("\tLoad graph from <load_dir>, then run the simulation.\n");
        return 0;
    }

    boost::filesystem::path dir1(result_dir);
    bool created_result_dir = boost::filesystem::create_directories(dir1);
    boost::filesystem::path dir2(save_dir);
    bool created_save_dir = boost::filesystem::create_directories(dir2);

    if (created_result_dir) {
        mio::log_info("Result directory '{:s}' was created.", dir1.string());
    }
    if (created_save_dir) {
        mio::log_info("Save directory '{:s}' was created.", dir2.string());
    }
    printf("Saving results to \"%s\".\n", result_dir.c_str());

    auto result = run(mode, data_dir, save_dir, result_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    return 0;
}
