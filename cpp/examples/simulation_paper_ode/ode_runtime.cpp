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
#include "memilio/config.h"
#include "models/ode_secir/infection_state.h"
#include "models/ode_secir/model.h"
#include "ode_secir/parameters_io.h"
#include "ode_secir/parameter_space.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/table_printer.h"
#include "memilio/utils/base_dir.h"
#include "memilio/io/cli.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/filesystem.hpp"
#include <iostream>
#include <vector>
#include <omp.h>


using Model    = mio::osecir::Model<ScalarType>;

mio::UncertainValue<ScalarType> uncertain(ScalarType v)
{
    const double var = .1;
    return mio::UncertainValue<ScalarType>(v, mio::ParameterDistributionUniform(v * (1 - var), v * (1 + var)));
}

namespace params
{
const ScalarType total_pop = 100000;
// Transition probabilities
const ScalarType symptomsPerInfectedNoSymptoms = 0.75;
const ScalarType severePerInfectedSymptoms     = 0.01;
const ScalarType criticalPerInfectedSevere     = 0.1;
const ScalarType deathsPerCritical             = 0.05;

// Mean stay times
const ScalarType timeExposed            = 4.5;
const ScalarType timeInfectedNoSymptoms = 3.;
const ScalarType timeInfectedSymptoms   = 7.;
const ScalarType timeInfectedSevere     = 15.;
const ScalarType timeInfectedCritical   = 16.;

const ScalarType transmissionProbabilityOnContact = 0.5;

// Mobility
const int band_radius_mobility_matrix = 50;
const ScalarType factorMobilePopulation = 0.1;

// Simulation parameters
ScalarType t0   = 0.;
ScalarType tmax = 30;
ScalarType dt   = 0.1;

} // namespace params

Model initialize_osecir()
{
    Model model(1);

    model.populations.set_total(params::total_pop);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = uncertain(0.01) * params::total_pop;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                params::total_pop);

    model.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[mio::AgeGroup(0)] = uncertain(params::timeExposed);
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(params::timeInfectedNoSymptoms);
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(params::timeInfectedSymptoms);
    model.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[mio::AgeGroup(0)] = uncertain(params::timeInfectedSevere);
    model.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(params::timeInfectedCritical);
    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(params::transmissionProbabilityOnContact);
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(1. - params::symptomsPerInfectedNoSymptoms);
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(params::severePerInfectedSymptoms);
    model.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[mio::AgeGroup(0)] =
        uncertain(params::criticalPerInfectedSevere);
    model.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[mio::AgeGroup(0)] = uncertain(params::deathsPerCritical);
    model.apply_constraints();

    mio::ContactMatrixGroup<ScalarType>& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                       = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, 1.0));

    model.apply_constraints();
    
    return model;
}

mio::Graph<Model, mio::MobilityParameters<ScalarType>> get_graph(size_t num_regions)
{
    mio::Graph<Model, mio::MobilityParameters<ScalarType>> params_graph;
    auto model = initialize_osecir();

    for (size_t region = 0; region < num_regions; region++)
    {
        params_graph.add_node((int)region, model);
    }

    // check for higher requested edges than number of nodes, then make it fully connected
    if ((params::band_radius_mobility_matrix * 2) >= num_regions){
        for (size_t region_out = 0; region_out < num_regions; region_out++)
        {   
            ScalarType mobilityPerEdge = params::factorMobilePopulation / (num_regions - 1);
            for (size_t region_in = 0; region_in < num_regions; region_in++)
            {
                if(region_out != region_in)
                    params_graph.add_edge(region_out, region_in, Eigen::VectorX<ScalarType>::Constant((size_t)mio::osecir::InfectionState::Count, mobilityPerEdge));
            }
        }
    }
    else {
        ScalarType mobilityPerEdge = params::factorMobilePopulation / (params::band_radius_mobility_matrix * 2);
        for (size_t region_out = 0; region_out < num_regions; region_out++)
        {   
            for (size_t idx_in = 1; idx_in < params::band_radius_mobility_matrix + 1; idx_in++)
            {
                size_t left  = (region_out + (int)num_regions - idx_in) % (int)num_regions;
                size_t right = (region_out + idx_in) % (int)num_regions;

                params_graph.add_edge(region_out, left, Eigen::VectorX<ScalarType>::Constant((size_t)mio::osecir::InfectionState::Count, mobilityPerEdge));
                params_graph.add_edge(region_out, right, Eigen::VectorX<ScalarType>::Constant((size_t)mio::osecir::InfectionState::Count, mobilityPerEdge));
            }
        }
    }
    return params_graph;
}

auto sample_graph(mio::Graph<Model, mio::MobilityParameters<ScalarType>> params_graph)
{
    auto copy = params_graph;
    return mio::make_sampled_graph_simulation<ScalarType, mio::osecir::Simulation<ScalarType>>(
        mio::osecir::draw_sample(copy), params::t0, params::dt, 0.5);
}

mio::IOResult<void> simulate_ode_secir(std::string result_dir, size_t num_runs, size_t num_warm_up_runs, size_t num_regions)
{
    // std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_result(
    //     num_runs,
    //     std::vector<mio::TimeSeries<ScalarType>>(1, mio::TimeSeries<ScalarType>(int(mio::osecir::InfectionState::Count))));

    std::vector<ScalarType> init_time(num_runs);
    std::vector<ScalarType> sim_time(num_runs);

    auto params_graph = get_graph(num_regions);

    mio::set_log_level(mio::LogLevel::off);
    for (size_t i = 0; i < num_warm_up_runs; i++) {
        auto sim = sample_graph(params_graph);
        sim.advance(params::tmax);
    }

    for (size_t run = 0; run < num_runs; run++) {
        mio::timing::BasicTimer timer;
        // Initialization
        timer.start();
        auto sim = sample_graph(params_graph);
        timer.stop();
        init_time[run] = timer.get_elapsed_time();
        timer.reset();
        //Simulation
        timer.start();
        sim.advance(params::tmax);
        timer.stop();
        sim_time[run] = timer.get_elapsed_time();
        // ensemble_result[run][0] = mio::interpolate_simulation_result(sim.get_result());
    }

    // Save ensemble results
    // BOOST_OUTCOME_TRY(save_results(ensemble_result, result_dir, "osecir"));
    mio::TimeSeries<ScalarType> init_time_ts(1);
    mio::TimeSeries<ScalarType> sim_time_ts(1);
    for (size_t i = 0; i < num_runs; i++) {
        init_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, init_time[i]));
        sim_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, sim_time[i]));
    }
    BOOST_OUTCOME_TRY(init_time_ts.export_csv(mio::path_join(result_dir, "osecir_" + std::to_string(num_regions) + "regions_init_time.csv")));
    BOOST_OUTCOME_TRY(sim_time_ts.export_csv(mio::path_join(result_dir,  "osecir_" + std::to_string(num_regions) + "regions_sim_time.csv")));

    return mio::success();
}

// void simulate_ode_secir(size_t num_runs, size_t num_warm_up_runs, size_t num_regions)
// {
//     using namespace params;

//     auto model = initialize_osecir();

//     // Perform simulation several times and measure run times.
//     // Warm up runs.
//     mio::set_log_level(mio::LogLevel::off);
//     for (size_t i = 0; i < num_warm_up_runs; i++) {
//         mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
//     }
//     // Simulate one time to track the number of steps.
//     auto result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
//     std::cout << "\"Steps\": " << result.get_num_time_points() << "," << std::endl;

//     // Runs with timing.
//     ScalarType total = 0;
//     for (size_t i = 0; i < num_runs; i++) {
//         total -= omp_get_wtime();
//         mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
//         total += omp_get_wtime();
//     }
//     std::cout << "\"Time\": " << total / num_runs << "\n}," << std::endl;
//     mio::set_log_level(mio::LogLevel::warn);
// }

int main(int argc, char** argv)
{
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                          .add<"ResultDirectory">(mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper_ode/results_runtime"))
                          .add<"NumberRuns">(100, {.alias = "nRun"})
                          .add<"NumberWarmupRuns">(10, {.alias = "nWURun"})
                          .add<"NumberRegions">(1000, {.alias = "nRegion"})
                          .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"ResultDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }

    boost::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    boost::filesystem::create_directories(res_dir);

    auto result = simulate_ode_secir(cli_parameters.get<"ResultDirectory">(), cli_parameters.get<"NumberRuns">(), cli_parameters.get<"NumberWarmupRuns">(), cli_parameters.get<"NumberRegions">());
    return 0;
}
