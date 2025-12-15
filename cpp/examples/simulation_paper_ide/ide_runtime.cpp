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
#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/simulation.h"

#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/utils/base_dir.h"
#include "memilio/utils/logging.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/table_printer.h"
#include "memilio/io/cli.h"

ScalarType uncertain(ScalarType v)
{
    const double var = .1;
    return mio::ParameterDistributionUniform(v * (1 - var), v * (1 + var)).get_sample(mio::thread_local_rng());
}

namespace params
{

// Mobility
const ScalarType factorMobilePopulation = 0.1;

// Simulation parameters
ScalarType t0   = 0.;
ScalarType tmax = 60;
ScalarType dt   = 0.1;

}

using Model    = mio::isecir::Model;

Model initialize_isecir()
{
    // Define number of age groups.
    size_t num_agegroups = 1;

    // Define initial values for the total population and number of deaths.
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), uncertain(10000.));
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 13.);

    // Create TimeSeries with num_transitions * num_agegroups elements where initial transitions needed for simulation
    // will be stored. We require values for the transitions for a sufficient number of time points before the start of
    // the simulation to initialize our model.
    size_t num_transitions = (size_t)mio::isecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> transitions_init(num_transitions * num_agegroups);

    // Define vector of transitions that will be added to the time points of the TimeSeries transitions_init.
    Eigen::VectorX<ScalarType> vec_init(num_transitions * num_agegroups);
    vec_init[(size_t)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = uncertain(25.);
    vec_init[(size_t)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = uncertain(15.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = uncertain(8.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = uncertain(4.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = uncertain(1.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = uncertain(4.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = uncertain(1.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = uncertain(1.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = uncertain(1.);
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = uncertain(1.);

    // define exponential function with memilio-simulation ide_convergence_rate.cpp anna2024
    // Multiply vec_init with dt_ide_solver so that within a time interval of length 1, always the a
    vec_init = vec_init * params::dt;
    transitions_init.add_time_point(-10, vec_init);
    while (transitions_init.get_last_time() < params::t0 - params::dt / 2.) {
        transitions_init.add_time_point(transitions_init.get_last_time() + params::dt, vec_init);
    }

    Model model(std::move(transitions_init), total_population_init, deaths_init, num_agegroups);
    model.check_constraints(params::dt);

    return model;
}

mio::Graph<mio::SimulationNode<ScalarType, mio::isecir::Simulation>, mio::MobilityEdge<ScalarType>> get_graph(size_t num_regions)
{
    mio::Graph<mio::SimulationNode<ScalarType, mio::isecir::Simulation>, mio::MobilityEdge<ScalarType>> sim_graph;

    for (size_t region = 0; region < num_regions; region++)
    {
        Model model = initialize_isecir();
        sim_graph.add_node((int)region, model, params::dt);
    }

    // ScalarType num_nodes = sim_graph.nodes().size();
    // for (size_t region_out = 0; region_out < num_regions; region_out++)
    // {   
    //     ScalarType mobilityPerEdge = params::factorMobilePopulation / num_nodes;
    //     for (size_t region_in = 0; region_in < num_regions; region_in++)
    //     {
    //         if (region_out != region_in)
    //         {
    //             sim_graph.add_edge(region_out, region_in, Eigen::VectorX<ScalarType>::Constant((size_t)mio::osecir::InfectionState::Count, mobilityPerEdge));
    //         }
            
    //     }
    // }
    return sim_graph;
}


mio::IOResult<void> simulate(std::string result_dir, size_t num_runs, size_t num_warm_up_runs, size_t num_regions)
{
    std::vector<ScalarType> init_time(num_runs);
    std::vector<ScalarType> sim_time(num_runs);

    mio::set_log_level(mio::LogLevel::off);
    for (size_t i = 0; i < num_warm_up_runs; i++) {
        auto sim_graph = get_graph(num_regions);
        auto sim = mio::make_no_mobility_sim(params::t0, std::move(sim_graph));
        sim.advance(params::tmax);
    }

    for (size_t run = 0; run < num_runs; run++) {
        mio::timing::BasicTimer timer;
        // Initialization
        timer.start();
        auto sim_graph = get_graph(num_regions);
        auto sim = mio::make_no_mobility_sim(params::t0, std::move(sim_graph));
        timer.stop();
        init_time[run] = timer.get_elapsed_time();
        timer.reset();
        //Simulation
        timer.start();
        sim.advance(params::tmax);
        timer.stop();
        sim_time[run]           = timer.get_elapsed_time();
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
    BOOST_OUTCOME_TRY(init_time_ts.export_csv(mio::path_join(result_dir, "isecir_" + std::to_string(num_regions) + "regions_init_time.csv")));
    BOOST_OUTCOME_TRY(sim_time_ts.export_csv(mio::path_join(result_dir,  "isecir_" + std::to_string(num_regions) + "regions_sim_time.csv")));

    return mio::success();
}

int main(int argc, char** argv)
{
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                          .add<"ResultDirectory">(mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper_ide/results_runtime"))
                          .add<"NumberRuns">(100, {.alias = "nRun"})
                          .add<"NumberWarmupRuns">(10, {.alias = "nWURun"})
                          .add<"NumberRegions">(10, {.alias = "nRegion"})
                          .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"ResultDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }

    boost::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    boost::filesystem::create_directories(res_dir);

    auto result = simulate(cli_parameters.get<"ResultDirectory">(), cli_parameters.get<"NumberRuns">(), cli_parameters.get<"NumberWarmupRuns">(), cli_parameters.get<"NumberRegions">());
    return 0;
}