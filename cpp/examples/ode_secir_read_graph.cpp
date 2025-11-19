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
#include "memilio/compartments/parameter_studies.h"
#include "memilio/config.h"
#include "memilio/io/cli.h"
#include "memilio/io/mobility_io.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/base_dir.h"

#include "memilio/utils/stl_util.h"
#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters_io.h"

#include <iostream>

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::critical);

    auto parameters =
        mio::cli::ParameterSetBuilder()
            .add<"MobilityFile">(
                mio::path_join(mio::base_dir(), "data", "Germany", "mobility", "commuter_mobility_2022.txt"),
                {.description = "Create the mobility file with MEmilio Epidata's getCommuterMobility.py file."})
            .build();

    auto result = mio::command_line_interface(argv[0], argc, argv, parameters, {"MobilityFile"});
    if (!result) {
        std::cout << result.error().message();
        return result.error().code().value();
    }

    const auto t0   = 0.;
    const auto tmax = 10.;

    ScalarType cont_freq = 10; // see Polymod study

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<ScalarType> model(1);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    ScalarType fact         = 1.0 / (ScalarType)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity<ScalarType>>(std::numeric_limits<ScalarType>::max());
    params.set<mio::osecir::StartDay<ScalarType>>(0);
    params.set<mio::osecir::Seasonality<ScalarType>>(0);

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<ScalarType>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 2.0;
        params.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 6.;
        params.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[i]     = 12;
        params.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[i]   = 8;

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

        params.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<ScalarType>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<ScalarType>>()[i]                = 0.3;
    }
    // The function apply_constraints() ensures that all parameters are within their defined bounds.
    // Note that negative values are set to zero instead of stopping the simulation.
    params.apply_constraints();

    mio::ContactMatrixGroup<ScalarType>& contact_matrix = params.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                                   = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    std::cout << "Reading Mobility File..." << std::flush;
    auto read_mobility_result = mio::read_mobility_plain(parameters.get<"MobilityFile">());
    if (!read_mobility_result) {
        std::cout << "\n" << read_mobility_result.error().formatted_message() << "\n";
        std::cout << "Create the mobility file with MEmilio Epidata's getCommuterMobility.py file.\n";
        return 0;
    }
    auto& commuter_mobility = read_mobility_result.value();
    std::cout << "Done" << std::endl;

    std::cout << "Intializing Graph..." << std::flush;
    mio::Graph<mio::osecir::Model<ScalarType>, mio::MobilityParameters<ScalarType>> graph;
    for (int node = 0; node < commuter_mobility.rows(); node++) {
        graph.add_node(node, model);
    }
    for (int row = 0; row < commuter_mobility.rows(); row++) {
        for (int col = 0; col < commuter_mobility.cols(); col++) {
            graph.add_edge(row, col,
                           Eigen::VectorX<ScalarType>::Constant(
                               10 * (size_t)nb_groups,
                               commuter_mobility(row, col) / graph.nodes()[row].property.populations.get_total()));
        }
    }
    std::cout << "Done" << std::endl;

    std::cout << "Writing Json Files..." << std::flush;
    auto write_status = mio::write_graph(graph, "graph_parameters");
    if (!write_status) {
        std::cout << "\n" << write_status.error().formatted_message();
        return 0;
    }
    std::cout << "Done" << std::endl;

    std::cout << "Reading Json Files..." << std::flush;
    auto graph_read_result = mio::read_graph<ScalarType, mio::osecir::Model<ScalarType>>("graph_parameters");

    if (!graph_read_result) {
        std::cout << "\n" << graph_read_result.error().formatted_message();
        return 0;
    }
    std::cout << "Done" << std::endl;
    auto& graph_read = graph_read_result.value();

    std::cout << "Running Simulations..." << std::flush;
    mio::ParameterStudy study(graph_read, t0, tmax, 0.5, 2);
    study.run_serial([](auto&& g, auto t0_, auto dt_, auto) {
        auto copy = g;
        return mio::make_sampled_graph_simulation<double, mio::osecir::Simulation<ScalarType>>(draw_sample(copy), t0_,
                                                                                               dt_, dt_);
    });
    std::cout << "Done" << std::endl;

    return 0;
}
