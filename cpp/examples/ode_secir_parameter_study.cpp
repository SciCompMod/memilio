/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "ode_secir/parameters_io.h"
#include "ode_secir/parameter_space.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/io/result_io.h"

/**
 * @brief creates xml file with a single run parameter study with std 0 (used to save parameters of individual runs)
 * @param filename Name of file
 * @param params Secir parameters used during run
 * @param t0 starting point of simulation
 * @param tmax end point of simulation
 */
mio::IOResult<void>
write_single_run_result(const int run,
                        const mio::Graph<mio::SimulationNode<mio::osecir::Simulation<>>, mio::MigrationEdge>& graph)
{
    std::string abs_path;
    BOOST_OUTCOME_TRY(created, mio::create_directory("results", abs_path));

    if (run == 0) {
        std::cout << "Results are stored in " << abs_path << '\n';
        if (!created) {
            std::cout << "Directory already exists, files from previous runs will be overwritten." << '\n';
        }
    }

    //write sampled parameters for this run
    //omit edges to save space as they are not sampled
    int inode = 0;
    for (auto&& node : graph.nodes()) {
        BOOST_OUTCOME_TRY(js_node_model, serialize_json(node.property.get_result(), mio::IOF_OmitDistributions));
        Json::Value js_node(Json::objectValue);
        js_node["NodeId"]  = node.id;
        js_node["Model"]   = js_node_model;
        auto node_filename = mio::path_join(abs_path, "Parameters_run" + std::to_string(run) + "_node" +
                                                          std::to_string(inode++) + ".json");
        BOOST_OUTCOME_TRY(mio::write_json(node_filename, js_node));
    }

    //write results for this run
    std::vector<mio::TimeSeries<double>> all_results;
    std::vector<int> ids;

    ids.reserve(graph.nodes().size());
    all_results.reserve(graph.nodes().size());
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(all_results), [](auto& node) {
        return node.property.get_result();
    });
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(ids), [](auto& node) {
        return node.id;
    });
    auto num_groups = (int)(size_t)graph.nodes()[0].property.get_simulation().get_model().parameters.get_num_groups();
    BOOST_OUTCOME_TRY(mio::save_result(all_results, ids, num_groups,
                                       mio::path_join(abs_path, ("Results_run" + std::to_string(run) + ".h5"))));

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;

    double cont_freq = 10; // see Polymod study

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::osecir::Model model(1);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<mio::osecir::StartDay>(0);
    params.set<mio::osecir::Seasonality>(0);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 6.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 12;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto write_parameters_status = mio::write_json("parameters.json", model);
    if (!write_parameters_status) {
        std::cout << "Error writing parameters: " << write_parameters_status.error().formatted_message();
        return -1;
    }

    //create study
    auto num_runs = size_t(1);
    mio::ParameterStudy<mio::osecir::Simulation<>> parameter_study(model, t0, tmax, num_runs);

    //run study
    auto sample_graph = [](auto&& graph) {
        return mio::osecir::draw_sample(graph);
    };
    auto handle_result = [](auto&& graph, auto&& run) {
        auto write_result_status = write_single_run_result(run, graph);
        if (!write_result_status) {
            std::cout << "Error writing result: " << write_result_status.error().formatted_message();
        }       
        return 0; //Result handler must return something, but only meaningful when using MPI.
    };
    parameter_study.run(sample_graph, handle_result);

    return 0;
}
