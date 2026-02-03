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
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/table_printer.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/utils/base_dir.h"
#include "memilio/io/cli.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <iostream>
#include <vector>
#include <omp.h>

mio::UncertainValue<ScalarType> uncertain(ScalarType v)
{
    const double var = .1;
    return mio::UncertainValue<ScalarType>(v, mio::ParameterDistributionUniform(v * (1 - var), v * (1 + var)));
}

namespace params
{
// num_subcompartments is used as a template argument and has to be a constexpr.
constexpr size_t num_subcompartments = 5;

// Define (non-age-resolved) parameters.
const ScalarType seasonality                    = 0.;
const ScalarType relativeTransmissionNoSymptoms = 1.;
const ScalarType riskOfInfectionFromSymptomatic = 0.3;
// const ScalarType total_population               = 83155031.0;

const ScalarType timeExposed                      = 3.335;
const ScalarType timeInfectedNoSymptoms           = 2.58916;
const ScalarType timeInfectedSymptoms             = 6.94547;
const ScalarType timeInfectedSevere               = 7.28196;
const ScalarType timeInfectedCritical             = 13.066;
const ScalarType transmissionProbabilityOnContact = 0.07333;
const ScalarType recoveredPerInfectedNoSymptoms   = 0.206901;
const ScalarType severePerInfectedSymptoms        = 0.07864;
const ScalarType criticalPerSevere                = 0.17318;
const ScalarType deathsPerCritical                = 0.21718;

// Vector is a "random vector" taken from another example. Just need some realistic values.
std::vector<ScalarType> init = {3966564.2110, 664.2367, 545.5523, 1050.3946, 5.6045, 0.5844, 307.4165, 0.};

// Mobility
const int band_radius_mobility_matrix = 50;
const ScalarType factorMobilePopulation = 0.1;

// Simulation parameters
ScalarType t0   = 0.;
ScalarType tmax = 30;
ScalarType dt   = 0.1;

} // namespace params

using InfState = mio::lsecir::InfectionState;
using LctState = mio::LctInfectionState<ScalarType, InfState, 1, params::num_subcompartments, params::num_subcompartments, params::num_subcompartments,
                                        params::num_subcompartments, params::num_subcompartments, 1, 1>;
using Model    = mio::lsecir::Model<ScalarType, LctState>;

Model initialize_lsecir()
{
    std::cout << "\"Subcompartments\": " << params::num_subcompartments << std::endl;
    // Initialize (non-age-resolved) LCT model.
    Model model;

    // Set parameters.
    model.parameters.get<mio::lsecir::TimeExposed<ScalarType>>()[0]                      = uncertain(params::timeExposed);
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0]           = uncertain(params::timeInfectedNoSymptoms);
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]             = uncertain(params::timeInfectedSymptoms);
    model.parameters.get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]               = uncertain(params::timeInfectedSevere);
    model.parameters.get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]             = uncertain(params::timeInfectedCritical);
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = uncertain(params::transmissionProbabilityOnContact);

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = uncertain(params::relativeTransmissionNoSymptoms);
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = uncertain(params::riskOfInfectionFromSymptomatic);
    model.parameters.get<mio::lsecir::Seasonality<ScalarType>>()                       = params::seasonality;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = uncertain(params::recoveredPerInfectedNoSymptoms);
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0]      = uncertain(params::severePerInfectedSymptoms);
    model.parameters.get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]              = uncertain(params::criticalPerSevere);
    model.parameters.get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]              = uncertain(params::deathsPerCritical);
    // Realistic average number of contacts.
    mio::ContactMatrixGroup<ScalarType> contact_matrix               = mio::ContactMatrixGroup<ScalarType>(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, 7.69129));
    model.parameters.get<mio::lsecir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix<ScalarType>(contact_matrix);

    // Use init as a basis to define appropriate initial values.
    // Compartment values are distributed equally to subcompartments.
    model.populations[0]                   = uncertain(params::init[0]); // Susceptible.
    model.populations[LctState::Count - 2] = params::init[6]; // Recovered.
    model.populations[LctState::Count - 1] = params::init[7]; // Dead.
    for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
        for (size_t subcomp = 0; subcomp < params::num_subcompartments; subcomp++) {
            model.populations[(i - 1) * params::num_subcompartments + 1 + subcomp] = uncertain(params::init[i] / (ScalarType)params::num_subcompartments);
        }
    }
    
    return model;
}

mio::Graph<Model, mio::MobilityParameters<ScalarType>> get_graph(size_t num_regions)
{
    mio::Graph<Model, mio::MobilityParameters<ScalarType>> params_graph;
    auto model = initialize_lsecir();

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
                    params_graph.add_edge(region_out, region_in, Eigen::VectorX<ScalarType>::Constant((size_t)LctState::Count, mobilityPerEdge));
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

                params_graph.add_edge(region_out, left, Eigen::VectorX<ScalarType>::Constant((size_t)LctState::Count, mobilityPerEdge));
                params_graph.add_edge(region_out, right, Eigen::VectorX<ScalarType>::Constant((size_t)LctState::Count, mobilityPerEdge));
            }
        }
    }
    return params_graph;
}

void draw_sample_infection(Model& model)
{
    //not age dependent
    model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0].draw_sample();
    model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0].draw_sample();
    model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0].draw_sample();
    model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0].draw_sample();
    model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0].draw_sample();

    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {
        //not age dependent
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[i] =
            model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0];
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[i] =
            model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0];
        model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[i] =
            model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0];
        model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[i] =
            model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0];

        //age dependent
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[i].draw_sample();
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[i].draw_sample();
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[i].draw_sample();

        model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[i].draw_sample();
        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i].draw_sample();
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[i].draw_sample();
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[i].draw_sample();
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[i].draw_sample();
    }
}

void draw_sample_demographics(Model& model)
{

    model.populations[LctState::Count - 2].draw_sample(); // Recovered.
    for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
        for (size_t subcomp = 0; subcomp < params::num_subcompartments; subcomp++) {
            model.populations[(i - 1) * params::num_subcompartments + 1 + subcomp].draw_sample();
        }
    }
}

mio::Graph<Model, mio::MobilityParameters<ScalarType>> draw_sample(mio::Graph<Model, mio::MobilityParameters<ScalarType>>& graph)
{
    mio::Graph<Model, mio::MobilityParameters<ScalarType>> sampled_graph;

    //sample global parameters
    auto& shared_params_model = graph.nodes()[0].property;
    draw_sample_infection(shared_params_model);
    auto& shared_contacts = shared_params_model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>();
    shared_contacts.draw_sample_dampings();

    for (auto& params_node : graph.nodes()) {
        auto& node_model = params_node.property;

        //sample local parameters
        draw_sample_demographics(params_node.property);

        node_model.parameters = shared_params_model.parameters;

        node_model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>().make_matrix();
        node_model.apply_constraints();

        sampled_graph.add_node(params_node.id, node_model);
    }

    for (auto& edge : graph.edges()) {
        auto edge_params = edge.property;
        mio::apply_dampings(edge_params.get_coefficients(), shared_contacts.get_dampings(), [&edge_params](auto& v) {
            return mio::make_mobility_damping_vector(edge_params.get_coefficients().get_shape(), v);
        });
        sampled_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge_params);
    }

    return sampled_graph;
}

auto sample_graph(mio::Graph<Model, mio::MobilityParameters<ScalarType>> params_graph)
{
    auto copy = params_graph;
    return mio::make_sampled_graph_simulation<ScalarType, mio::Simulation<ScalarType, Model>>(
        draw_sample(copy), params::t0, params::dt, 0.5);
}

mio::IOResult<void> simulate(std::string result_dir, size_t num_runs, size_t num_warm_up_runs, size_t num_regions)
{
    // // Set integrator of fifth order.
    // auto integrator =
    //     std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // // Choose dt_min = dt_max to get a fixed step size.
    // integrator->set_dt_min(dt);
    // integrator->set_dt_max(dt);

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
    BOOST_OUTCOME_TRY(init_time_ts.export_csv(mio::path_join(result_dir, "lsecir_" + std::to_string(num_regions) + "regions_init_time.csv")));
    BOOST_OUTCOME_TRY(sim_time_ts.export_csv(mio::path_join(result_dir,  "lsecir_" + std::to_string(num_regions) + "regions_sim_time.csv")));

    return mio::success();
}

int main(int argc, char** argv)
{
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                          .add<"ResultDirectory">(mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper_lct/results_runtime"))
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
