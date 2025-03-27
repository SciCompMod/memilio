/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/io/mobility_io.h"
#include "memilio/compartments/parameter_studies.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters_io.h"
#include <data_dir.h>
#include <iostream>

std::string setup(int argc, char** argv, const std::string data_dir)
{
    if (argc == 2) {
        std::cout << "Using file " << argv[1] << " in data/Germany/mobility." << std::endl;
        return mio::path_join(data_dir, "Germany", "mobility", (std::string)argv[1]);
    }
    else {
        if (argc > 2) {
            mio::log_error("Too many arguments given.");
        }
        else {
            mio::log_warning("No arguments given.");
        }
        auto mobility_file = "commuter_mobility_2022.txt";
        std::cout << "Using file " << mobility_file << " in data/Germany/mobility." << std::endl;
        std::cout << "Usage: read_graph MOBILITY_FILE"
                  << "\n\n";
        std::cout
            << "This example performs a simulation based on mobility data from the German Federal Employment Agency."
            << std::endl;
        return mio::path_join(data_dir, "Germany", "mobility", mobility_file);
    }
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::critical);
    std::string data_dir = DATA_DIR;
    std::string filename = setup(argc, argv, data_dir);

    const auto t0   = 0.;
    const auto tmax = 10.;
    using FP        = double;

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<FP> model(1);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact             = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity<FP>>(std::numeric_limits<double>::max());
    params.set<mio::osecir::StartDay>(0);
    params.set<mio::osecir::Seasonality<FP>>(0);

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.0;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 6.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 12;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<FP>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<FP>>()[i]                = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<FP>>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    std::cout << "Reading Mobility File..." << std::flush;
    auto read_mobility_result = mio::read_mobility_plain(filename);
    if (!read_mobility_result) {
        std::cout << read_mobility_result.error().formatted_message() << '\n';
        std::cout << "Create the mobility file with MEmilio Epidata's getCommuterMobility.py file." << '\n';
        return 0;
    }
    auto& commuter_mobility = read_mobility_result.value();
    std::cout << "Done" << std::endl;

    std::cout << "Intializing Graph..." << std::flush;
    mio::Graph<mio::osecir::Model<FP>, mio::MobilityParameters<FP>> graph;
    for (int node = 0; node < commuter_mobility.rows(); node++) {
        graph.add_node(node, model);
    }
    for (int row = 0; row < commuter_mobility.rows(); row++) {
        for (int col = 0; col < commuter_mobility.cols(); col++) {
            graph.add_edge(row, col,
                           Eigen::VectorXd::Constant(10 * (size_t)nb_groups,
                                                     commuter_mobility(row, col) /
                                                         graph.nodes()[row].property.populations.get_total()));
        }
    }
    std::cout << "Done" << std::endl;

    std::cout << "Writing Json Files..." << std::flush;
    auto write_status = mio::write_graph(graph, "graph_parameters");
    if (!write_status) {
        std::cout << "Error: " << write_status.error().formatted_message();
    }
    std::cout << "Done" << std::endl;

    std::cout << "Reading Json Files..." << std::flush;
    auto graph_read_result = mio::read_graph<FP, mio::osecir::Model<FP>>("graph_parameters");

    if (!graph_read_result) {
        std::cout << "Error: " << graph_read_result.error().formatted_message();
    }
    std::cout << "Done" << std::endl;
    auto& graph_read = graph_read_result.value();

    std::cout << "Running Simulations..." << std::flush;
    auto study = mio::ParameterStudy<mio::osecir::Simulation<>>(graph_read, t0, tmax, 0.5, 2);
    std::cout << "Done" << std::endl;

    return 0;
}
