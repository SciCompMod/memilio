/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "secir/secir_parameters_io.h"
#include "secir/parameter_space.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/mobility/mobility.h"
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
                        const mio::Graph<mio::SimulationNode<mio::SecirSimulation<>>, mio::MigrationEdge>& graph)
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
    BOOST_OUTCOME_TRY(
        mio::save_result(all_results, ids, num_groups,
                    mio::path_join(abs_path, ("Results_run" + std::to_string(run) + ".h5"))));

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;

    double tinc    = 5.2, 
        tinf   = 6, 
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        tsevere = 12, 
        tcritical  = 8;

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    mio::SecirModel model(1);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    params.set<mio::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<mio::StartDay>(0);
    params.set<mio::Seasonality>(0);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::IncubationTime>()[i]         = tinc;
        params.get<mio::TimeInfectedSymptoms>()[i]     = tinf;
        params.get<mio::SerialInterval>()[i]         = tserint;
        params.get<mio::TimeInfectedSevere>()[i] = tsevere;
        params.get<mio::TimeInfectedCritical>()[i]          = tcritical;

        model.populations[{i, mio::InfectionState::Exposed}]      = fact * num_exp_t0;
        model.populations[{i, mio::InfectionState::Carrier}]      = fact * num_car_t0;
        model.populations[{i, mio::InfectionState::Infected}]     = fact * num_inf_t0;
        model.populations[{i, mio::InfectionState::Hospitalized}] = fact * num_hosp_t0;
        model.populations[{i, mio::InfectionState::ICU}]          = fact * num_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]    = fact * num_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]         = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<mio::RelativeTransmissionNoSymptoms>()[i]    = carr_infec;
        params.get<mio::AsymptomaticCasesPerInfectious>()[i]    = alpha;
        params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::HospitalizedCasesPerInfectious>()[i]  = rho;
        params.get<mio::ICUCasesPerHospitalized>()[i]         = theta;
        params.get<mio::DeathsPerICU>()[i]                    = delta;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto write_parameters_status = mio::write_json("parameters.json", model);
    if (!write_parameters_status) {
        std::cout << "Error writing parameters: " << write_parameters_status.error().formatted_message();
        return -1;
    }

    //create study
    auto num_runs = size_t(1);
    mio::ParameterStudy<mio::SecirSimulation<>> parameter_study(model, t0, tmax, num_runs);

    //run study
    int run           = 0;
    auto sample_graph = [](auto&& graph) {
        return mio::draw_sample(graph);
    };
    auto handle_result = [&run](auto&& graph) {
        auto write_result_status = write_single_run_result(run++, graph);
        if (!write_result_status) {
            std::cout << "Error writing result: " << write_result_status.error().formatted_message();
        }
    };
    parameter_study.run(sample_graph, handle_result);

    return 0;
}
