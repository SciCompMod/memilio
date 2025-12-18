/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Maximilian Betz
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
#include "models/ode_secir/infection_state.h"
#include "models/ode_secir/model.h"
#include "models/ode_secir/parameter_space.h"
#include "models/ode_secir/parameters_io.h"

#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/base_dir.h"
#include "memilio/io/cli.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"
#include "memilio/utils/miompi.h"

#include <omp.h>
#include <mpi.h>

using Model    = mio::osecir::Model<ScalarType>;

mio::UncertainValue<ScalarType> uncertain(ScalarType v)
{
    const ScalarType var = .1;
    return mio::UncertainValue<ScalarType>(v, mio::ParameterDistributionUniform(v * (1 - var), v * (1 + var)));
}

namespace params
{
    constexpr size_t num_groups       = 6;
    size_t num_regions                = 400;
    int num_processes = 1;

    mio::Date start_date(2021, 01, 01);
    const ScalarType seasonality              = 0.;
    ScalarType relativeTransmissionNoSymptoms = 1.;
    ScalarType riskOfInfectionFromSymptomatic = 0.3;
    const ScalarType transmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

    // Transition probabilities
    const ScalarType recoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
    const ScalarType severePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
    const ScalarType criticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
    const ScalarType deathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

    // Mean stay times
    const ScalarType timeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
    const ScalarType timeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
    const ScalarType timeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
    const ScalarType timeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
    const ScalarType timeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

    const ScalarType age_group_sizes[]        = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
    const ScalarType total_population         = 83155031.0;

    // Mobility
    const int band_radius_mobility_matrix = 50;
    const ScalarType factorMobilePopulation = 0.1;

    // Simulation parameters
    ScalarType t0   = 0.;
    ScalarType tmax = 30;
    ScalarType dt   = 0.1;

    enum class ContactLocation
    {
        Home = 0,
        School,
        Work,
        Other,
        Count
    };

    static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                            {ContactLocation::School, "school_pf_eig"},
                                                                            {ContactLocation::Work, "work"},
                                                                            {ContactLocation::Other, "other"}};

}

mio::IOResult<void> set_contact_matrices(std::string data_dir, mio::osecir::Parameters<ScalarType>& params)
{
    auto num_groups = size_t(params.get_num_groups());
    auto contact_matrices = mio::ContactMatrixGroup<ScalarType>(params::contact_locations.size(), num_groups);
    for (auto&& contact_location : params::contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(mio::path_join(data_dir, "Germany", "contacts", ("baseline_" + contact_location.second + ".txt"))));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline.cast<ScalarType>();
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Zero(num_groups, num_groups);
    }
    params.template get<mio::osecir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix<ScalarType>(contact_matrices);

    return mio::success();
}


mio::osecir::Parameters<ScalarType> initialize_osecir_params(std::string data_dir)
{
    mio::osecir::Parameters<ScalarType> parameters(params::num_groups);

    for (size_t group = 0; group < params::num_groups; group++) {
        parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::transmissionProbabilityOnContact[group]);
        parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::relativeTransmissionNoSymptoms);
        parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::riskOfInfectionFromSymptomatic);

        // Mean stay times
        parameters.get<mio::osecir::TimeExposed<ScalarType>>()[mio::AgeGroup(group)] = uncertain(params::timeExposed[group]);
        parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::timeInfectedNoSymptoms[group]);
        parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::timeInfectedSymptoms[group]);
        parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[mio::AgeGroup(group)] = uncertain(params::timeInfectedSevere[group]);
        parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::timeInfectedCritical[group]);

        // Transition probabilities
        parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::recoveredPerInfectedNoSymptoms[group]);
        parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::severePerInfectedSymptoms[group]);
        parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[mio::AgeGroup(group)] =
            uncertain(params::criticalPerSevere[group]);
        parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[mio::AgeGroup(group)] = uncertain(params::deathsPerCritical[group]);
    }
    
    parameters.get<mio::osecir::Seasonality<ScalarType>>() = params::seasonality;
    parameters.get<mio::osecir::StartDay<ScalarType>>()    = mio::get_day_in_year(params::start_date);

    auto result = set_contact_matrices(data_dir, parameters);
    
    return parameters;
}

mio::IOResult<mio::Graph<Model, mio::MobilityParameters<ScalarType>>> get_graph(std::string data_dir)
{
    auto model_params = initialize_osecir_params(data_dir);

    mio::Graph<mio::osecir::Model<ScalarType>, mio::MobilityParameters<ScalarType>> params_graph;
    const auto& read_function_nodes = mio::osecir::read_input_data_county<mio::osecir::Model<ScalarType>>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    auto end_date        = mio::offset_date_by_days(params::start_date, params::tmax);
    auto scaling_factor_infected = std::vector<ScalarType>(size_t(model_params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    auto mobile_compartments     = {mio::osecir::InfectionState::Susceptible,
                                    mio::osecir::InfectionState::Exposed,
                                    mio::osecir::InfectionState::InfectedNoSymptoms,
                                    mio::osecir::InfectionState::InfectedSymptoms,
                                    mio::osecir::InfectionState::Recovered};
    auto tnt_capacity_factor     = 0.;

    const auto& set_node_function =
        mio::set_nodes<ScalarType, mio::osecir::TestAndTraceCapacity<ScalarType>, mio::osecir::ContactPatterns<ScalarType>,
                       mio::osecir::Model<ScalarType>, mio::MobilityParameters<ScalarType>,
                       mio::osecir::Parameters<ScalarType>, decltype(read_function_nodes), decltype(node_id_function)>;
    const auto& set_edge_function =
        mio::set_edges<ScalarType, params::ContactLocation, mio::osecir::Model<ScalarType>, mio::MobilityParameters<ScalarType>,
                       mio::MobilityCoefficientGroup<ScalarType>, mio::osecir::InfectionState, decltype(read_function_edges)>;
    
    auto pydata_path =
        mio::path_join(data_dir, "Germany", "pydata");
    auto population_data_path =
        mio::path_join(pydata_path, "county_current_population.json");

    BOOST_OUTCOME_TRY(set_node_function(model_params, params::start_date, end_date, pydata_path, population_data_path, true,
                                 params_graph, read_function_nodes, node_id_function, scaling_factor_infected,
                                 scaling_factor_icu, tnt_capacity_factor, 0, false , true));
    BOOST_OUTCOME_TRY(set_edge_function(mio::path_join(data_dir, "Germany", "mobility", "commuter_mobility_2019.txt"), params_graph, mobile_compartments, 
                                        params::contact_locations.size(), read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.},
                                        std::vector<std::vector<size_t>>{}));

    return mio::success(params_graph);
}

mio::IOResult<void> simulate(std::string result_dir, std::string data_dir, size_t num_ensemble_runs)
{
    if (mio::mpi::is_root()) {
        std::cout << "Realistic scenario." << std::endl;
    }
    auto params_graph = get_graph(data_dir).value();

    mio::ParameterStudy parameter_study(params_graph, params::t0, params::tmax, params::dt, num_ensemble_runs);
    
    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : parameter_study.get_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }
    ScalarType total_time = 0;
    if (mio::mpi::is_root()) {
        total_time -= omp_get_wtime();
    }
    parameter_study.run(
        [](auto&& graph, ScalarType t0, ScalarType dt, size_t) {
            auto copy = graph;
            return mio::make_sampled_graph_simulation<ScalarType, mio::osecir::Simulation<ScalarType>>(
                mio::osecir::draw_sample(copy), t0, dt, 0.5);
        },
        [](auto&& results_sim, auto&& run_id){
            auto results_graph       = results_sim.get_graph();
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);

            // auto params = std::vector<mio::osecir::Model<ScalarType>>{};
            // params.reserve(results_graph.nodes().size());
            // std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
            //                              [](auto&& node) {
            //                    return node.property.get_simulation().get_model();
            //                });

            // auto edges = std::vector<mio::TimeSeries<ScalarType>>{};
            // edges.reserve(results_graph.edges().size());
            // std::transform(results_graph.edges().begin(), results_graph.edges().end(), std::back_inserter(edges),
            //                              [](auto&& edge) {
            //                    return edge.property.get_mobility_results();
            //                });
            mio::unused(interpolated_result, run_id);
            // return std::make_tuple(std::move(interpolated_result), std::move(params), std::move(edges));
        }
        );
    if (mio::mpi::is_root()) {
        total_time += omp_get_wtime();
    }

    if (mio::mpi::is_root()) {
        std::cout << "\"Processes\": " << params::num_processes << "," << std::endl;
        std::cout << "\"Num_ensemble_runs\": " << num_ensemble_runs << "," << std::endl;
        std::cout << "\"Time\": " << total_time << "\n}," << std::endl;
    }

    mio::unused(result_dir);
    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::err);
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                          .add<"ResultDirectory">(mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper_ode/results_ensemble"))
                          .add<"DataDirectory">(mio::path_join(mio::base_dir(), "data"))
                          .add<"NumberEnsembleRuns">(100, {.alias = "nRun"})
                          .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"ResultDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }

    mio::mpi::init();
    int size;
    MPI_Comm_size(mio::mpi::get_world(), &size);
    params::num_processes = size;

    boost::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    boost::filesystem::create_directories(res_dir);

    auto result = simulate(cli_parameters.get<"ResultDirectory">(), cli_parameters.get<"DataDirectory">(), cli_parameters.get<"NumberEnsembleRuns">());
    if (!result) {
        if (mio::mpi::is_root()) {
            printf("%s\n", result.error().formatted_message().c_str());
        }
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}