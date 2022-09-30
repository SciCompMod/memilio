/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "memilio/mobility/mobility.h"
#include "memilio/io/mobility_io.h"
#include "memilio/compartments/parameter_studies.h"
#include "secir/parameter_space.h"
#include "secir/secir_parameters_io.h"
#include <data_dir.h>
#include <iostream>

void print_usage()
{
    std::cout << "Usage: read_graph MIGRATION_FILE"
              << "\n\n";
    std::cout << "This example performs a simulation based on twitter "
                 "migration data."
              << std::endl;
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        print_usage();
        return -1;
    }

    std::string dir = DATA_DIR;

    auto filename   = mio::path_join(dir, "migration", (std::string)argv[1]);
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of migration, not integration

    double tinc    = 5.2, 
        tinf   = 6, 
        tserint    = 4.2,
        tsevere = 12, 
        tcritical  = 8;

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::SecirModel model(1);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact             = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<mio::StartDay>(0);
    params.set<mio::Seasonality>(0);

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::IncubationTime>()[i]         = tinc;
        params.get<mio::TimeInfectedSymptoms>()[i]     = tinf;
        params.get<mio::SerialInterval>()[i]         = tserint;
        params.get<mio::TimeInfectedSevere>()[i] = tsevere;
        params.get<mio::TimeInfectedCritical>()[i]          = tcritical;

        model.populations[{i, mio::InfectionState::Exposed}]      = fact * nb_exp_t0;
        model.populations[{i, mio::InfectionState::Carrier}]      = fact * nb_car_t0;
        model.populations[{i, mio::InfectionState::Infected}]     = fact * nb_inf_t0;
        model.populations[{i, mio::InfectionState::Hospitalized}] = fact * nb_hosp_t0;
        model.populations[{i, mio::InfectionState::ICU}]          = fact * nb_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]    = fact * nb_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]         = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<mio::RelativeCarrierInfectability>()[i]    = carr_infec;
        params.get<mio::AsymptomaticCasesPerInfectious>()[i]    = alpha;
        params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::HospitalizedCasesPerInfectious>()[i]  = rho;
        params.get<mio::ICUCasesPerHospitalized>()[i]         = theta;
        params.get<mio::DeathsPerICU>()[i]                    = delta;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    mio::set_params_distributions_normal(model, t0, tmax, 0.2);

    std::cout << "Reading Migration File..." << std::flush;
    auto read_mobility_result = mio::read_mobility_plain(filename);
    if (!read_mobility_result) {
        std::cout << read_mobility_result.error().formatted_message() << '\n';
        return -1;
    }
    auto& twitter_migration_2018 = read_mobility_result.value();
    std::cout << "Done" << std::endl;

    std::cout << "Intializing Graph..." << std::flush;
    mio::Graph<mio::SecirModel, mio::MigrationParameters> graph;
    for (int node = 0; node < twitter_migration_2018.rows(); node++) {
        graph.add_node(node, model);
    }
    for (int row = 0; row < twitter_migration_2018.rows(); row++) {
        for (int col = 0; col < twitter_migration_2018.cols(); col++) {
            graph.add_edge(row, col,
                           Eigen::VectorXd::Constant(8 * (size_t)nb_groups,
                                                     twitter_migration_2018(row, col) /
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
    auto graph_read_result = mio::read_graph<mio::SecirModel>("graph_parameters");

    if (!graph_read_result) {
        std::cout << "Error: " << graph_read_result.error().formatted_message();
    }
    std::cout << "Done" << std::endl;
    auto& graph_read = graph_read_result.value();

    std::cout << "Running Simulations..." << std::flush;
    auto study = mio::ParameterStudy<mio::SecirSimulation<>>(graph_read, t0, tmax, 0.5, 2);
    std::cout << "Done" << std::endl;

    return 0;
}
