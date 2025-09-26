/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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
#include "memilio/io/io.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/utils/parameter_set.h"
#include "models/d_abm/model.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "models/hybrid/temporal_hybrid_model.h"
#include "models/d_abm/simulation.h"
#include "models/d_abm/single_well.h"
#include "memilio/compartments/simulation.h"
#include "models/ode_secir/infection_state.h"
#include "models/ode_secir/model.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/epidemiology/adoption_rate.h"
#include "memilio/geography/regions.h"
#include "ode_secir/infection_state.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/table_printer.h"
#include "models/hybrid/conversion_functions.cpp"
#include <boost/outcome/try.hpp>
#include <compare>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <filesystem>

namespace params
{
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

// Simulation parameters
ScalarType t0   = 0.;
ScalarType tmax = 60;
ScalarType dt   = 0.1;

// Special dABM parameters
const ScalarType interaction_radius = 0.1;
const ScalarType noise              = 0.5;

// Special hybrid parameters
ScalarType dt_switch = 0.2;

} // namespace params

mio::dabm::Model<SingleWell<mio::osecir::InfectionState>> initialize_abm(int total_pop, double init_E)
{
    using Model = mio::dabm::Model<SingleWell<mio::osecir::InfectionState>>;

    std::vector<Model::Agent> agents(total_pop);
    //Random variables for initialization of agents' position
    auto& pos_sampler    = mio::UniformDistribution<double>::get_instance();
    auto& status_sampler = mio::DiscreteDistribution<size_t>::get_instance();
    for (auto& a : agents) {
        //Agents' positions are equally distributed in [-2, 2] x [-2, 2]
        a.position = Eigen::Vector2d{pos_sampler(mio::thread_local_rng(), -2., 2.),
                                     pos_sampler(mio::thread_local_rng(), -2., 2.)};
        //Sample agents status
        bool is_exposed = status_sampler(mio::thread_local_rng(),
                                         std::vector<double>{total_pop - init_E * total_pop, init_E * total_pop}) == 1;
        a.status        = is_exposed ? mio::osecir::InfectionState::Exposed : mio::osecir::InfectionState::Susceptible;
    }

    // Initialize adoption rates
    std::vector<mio::AdoptionRate<mio::osecir::InfectionState>> adoption_rates;
    //Second-order adoption rate (S->E)
    adoption_rates.push_back(
        {mio::osecir::InfectionState::Susceptible,
         mio::osecir::InfectionState::Exposed,
         mio::regions::Region(0),
         params::transmissionProbabilityOnContact,
         {{mio::osecir::InfectionState::InfectedNoSymptoms, 1}, {mio::osecir::InfectionState::InfectedSymptoms, 1}}});
    //First-order adoption rates
    //E->Ins
    adoption_rates.push_back({mio::osecir::InfectionState::Exposed,
                              mio::osecir::InfectionState::InfectedNoSymptoms,
                              mio::regions::Region(0),
                              1. / params::timeExposed,
                              {}});
    //Ins->Isy
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedNoSymptoms,
                              mio::osecir::InfectionState::InfectedSymptoms,
                              mio::regions::Region(0),
                              params::symptomsPerInfectedNoSymptoms / params::timeInfectedNoSymptoms,
                              {}});
    //Ins->R
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedNoSymptoms,
                              mio::osecir::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - params::symptomsPerInfectedNoSymptoms) / params::timeInfectedNoSymptoms,
                              {}});
    //Isy->Isev
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedSymptoms,
                              mio::osecir::InfectionState::InfectedSevere,
                              mio::regions::Region(0),
                              params::severePerInfectedSymptoms / params::timeInfectedSymptoms,
                              {}});
    //Isy->R
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedSymptoms,
                              mio::osecir::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - params::severePerInfectedSymptoms) / params::timeInfectedSymptoms,
                              {}});
    //Isev->Icri
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedSevere,
                              mio::osecir::InfectionState::InfectedCritical,
                              mio::regions::Region(0),
                              params::criticalPerInfectedSevere / params::timeInfectedSevere,
                              {}});
    //Isev->R
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedSevere,
                              mio::osecir::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - params::criticalPerInfectedSevere) / params::timeInfectedSevere,
                              {}});
    //Icri->R
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedCritical,
                              mio::osecir::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - params::deathsPerCritical) / params::timeInfectedCritical,
                              {}});
    //Icri->D
    adoption_rates.push_back({mio::osecir::InfectionState::InfectedCritical,
                              mio::osecir::InfectionState::Dead,
                              mio::regions::Region(0),
                              params::deathsPerCritical / params::timeInfectedCritical,
                              {}});

    // Create model
    Model model(agents, adoption_rates, params::interaction_radius, params::noise,
                {mio::osecir::InfectionState::InfectedSevere, mio::osecir::InfectionState::InfectedCritical,
                 mio::osecir::InfectionState::Dead});
    return model;
}

mio::IOResult<void> save_results(std::vector<std::vector<mio::TimeSeries<double>>>& ensemble_result,
                                 std::string result_dir, std::string model_name)
{
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.05)[0].export_csv(
        result_dir + model_name + "_p05.csv", {"S", "E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri", "R", "D"}));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.25)[0].export_csv(
        result_dir + model_name + "_p25.csv", {"S", "E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri", "R", "D"}));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.50)[0].export_csv(
        result_dir + model_name + "_p50.csv", {"S", "E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri", "R", "D"}));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.75)[0].export_csv(
        result_dir + model_name + "_p75.csv", {"S", "E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri", "R", "D"}));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.95)[0].export_csv(
        result_dir + model_name + "_p95.csv", {"S", "E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri", "R", "D"}));
    return mio::success();
}

mio::IOResult<void> simulate_dabm(std::string result_dir, size_t num_runs, int total_pop, double init_E)
{
    // As we only have one region, ensemble outputs are a vector of timeseries with size 1 for every run
    std::vector<std::vector<mio::TimeSeries<double>>> ensemble_result(
        num_runs,
        std::vector<mio::TimeSeries<double>>(1, mio::TimeSeries<double>(int(mio::osecir::InfectionState::Count))));

    std::vector<double> init_time(num_runs);
    std::vector<double> sim_time(num_runs);

#pragma omp parallel for
    for (size_t run = 0; run < num_runs; run++) {
        // Initialization
        mio::timing::BasicTimer timer;
        timer.start();
        auto model = initialize_abm(total_pop, init_E);
        auto sim   = mio::dabm::Simulation(model, params::t0, params::dt);
        timer.stop();
        init_time[run] = timer.get_elapsed_time();
        timer.reset();
        // Simulation
        timer.start();
        sim.advance(params::tmax);
        timer.stop();
        sim_time[run]           = timer.get_elapsed_time();
        ensemble_result[run][0] = mio::interpolate_simulation_result(sim.get_result());
    }
    // Save ensemble results
    BOOST_OUTCOME_TRY(save_results(ensemble_result, result_dir, "dabm"));
    // Convert init_time and sim_time to TimeSeries and save as csv
    mio::TimeSeries<double> init_time_ts(1);
    mio::TimeSeries<double> sim_time_ts(1);
    for (size_t i = 0; i < num_runs; i++) {
        init_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, init_time[i]));
        sim_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, sim_time[i]));
    }
    BOOST_OUTCOME_TRY(init_time_ts.export_csv(result_dir + "dabm_init_time.csv"));
    BOOST_OUTCOME_TRY(sim_time_ts.export_csv(result_dir + "dabm_sim_time.csv"));
    return mio::success();
}

mio::osecir::Model<double> initialize_osecir(int total_pop, double init_E)
{
    mio::osecir::Model<double> model(1);

    model.populations.set_total(total_pop);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = init_E * total_pop;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                total_pop);

    model.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[mio::AgeGroup(0)] = params::timeExposed;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        params::timeInfectedNoSymptoms;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        params::timeInfectedSymptoms;
    model.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[mio::AgeGroup(0)] = params::timeInfectedSevere;
    model.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[mio::AgeGroup(0)] =
        params::timeInfectedCritical;
    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(0)] =
        params::transmissionProbabilityOnContact;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        1. - params::symptomsPerInfectedNoSymptoms;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[mio::AgeGroup(0)] =
        params::severePerInfectedSymptoms;
    model.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[mio::AgeGroup(0)] =
        params::criticalPerInfectedSevere;
    model.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[mio::AgeGroup(0)] = params::deathsPerCritical;
    model.apply_constraints();
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.0));
    model.apply_constraints();
    return model;
}

mio::IOResult<void> simulate_ode_secir(std::string result_dir, size_t num_runs, int total_pop, double init_E)
{
    // As we don't have any spatial resolution, ensemble outputs are a vector of timeseries with size 1 for every run
    std::vector<std::vector<mio::TimeSeries<double>>> ensemble_result(
        num_runs,
        std::vector<mio::TimeSeries<double>>(1, mio::TimeSeries<double>(int(mio::osecir::InfectionState::Count))));

    std::vector<double> init_time(num_runs);
    std::vector<double> sim_time(num_runs);

#pragma omp parallel for
    for (size_t run = 0; run < num_runs; run++) {
        mio::timing::BasicTimer timer;
        // Initialization
        timer.start();
        auto model = initialize_osecir(total_pop, init_E);
        auto sim   = mio::Simulation(model, params::t0, params::dt);
        timer.stop();
        init_time[run] = timer.get_elapsed_time();
        timer.reset();
        //Simulation
        timer.start();
        sim.advance(params::tmax);
        timer.stop();
        sim_time[run]           = timer.get_elapsed_time();
        ensemble_result[run][0] = mio::interpolate_simulation_result(sim.get_result());
    }

    // Save ensemble results
    BOOST_OUTCOME_TRY(save_results(ensemble_result, result_dir, "osecir"));
    mio::TimeSeries<double> init_time_ts(1);
    mio::TimeSeries<double> sim_time_ts(1);
    for (size_t i = 0; i < num_runs; i++) {
        init_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, init_time[i]));
        sim_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, sim_time[i]));
    }
    BOOST_OUTCOME_TRY(init_time_ts.export_csv(result_dir + "osecir_init_time.csv"));
    BOOST_OUTCOME_TRY(sim_time_ts.export_csv(result_dir + "osecir_sim_time.csv"));
    return mio::success();
}

mio::IOResult<void> simulate_hybrid(std::string result_dir, size_t num_runs, const double switch_threshold,
                                    int total_pop, double init_E)
{

    std::vector<std::vector<mio::TimeSeries<double>>> ensemble_result(
        num_runs,
        std::vector<mio::TimeSeries<double>>(1, mio::TimeSeries<double>(int(mio::osecir::InfectionState::Count))));

    std::vector<double> init_time(num_runs);
    std::vector<double> sim_time(num_runs);

    size_t switch_runs = 0;

#pragma omp parallel for
    for (size_t run = 0; run < num_runs; run++) {
        mio::timing::BasicTimer timer;
        // Initialization
        timer.start();
        auto abm     = initialize_abm(total_pop, init_E);
        auto pbm     = initialize_osecir(total_pop, init_E);
        auto sim_abm = mio::dabm::Simulation(abm, params::t0, params::dt);
        auto sim_pbm = mio::Simulation(pbm, params::t0, params::dt);

        const auto result_fct_abm = [](const mio::dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& sim,
                                       double /*t*/) {
            return sim.get_result();
        };

        const auto result_fct_pbm = [](const mio::Simulation<double, mio::osecir::Model<double>>& sim, double /*t*/) {
            return sim.get_result();
        };

        //Define switching condition

        using HybridSim = mio::hybrid::TemporalHybridSimulation<decltype(sim_abm), decltype(sim_pbm),
                                                                mio::TimeSeries<double>, mio::TimeSeries<double>>;

        //Define switching condition
        const auto condition = [&switch_threshold, total_pop](const mio::TimeSeries<double>& result_abm,
                                                              const mio::TimeSeries<double>& /*result_pbm*/,
                                                              bool abm_used) {
            if (abm_used) {
                auto& last_value    = result_abm.get_last_value().eval();
                double num_infected = last_value[(int)mio::osecir::InfectionState::Exposed] +
                                      last_value[(int)mio::osecir::InfectionState::InfectedNoSymptoms] +
                                      last_value[(int)mio::osecir::InfectionState::InfectedNoSymptomsConfirmed] +
                                      last_value[(int)mio::osecir::InfectionState::InfectedSymptoms] +
                                      last_value[(int)mio::osecir::InfectionState::InfectedSymptomsConfirmed] +
                                      last_value[(int)mio::osecir::InfectionState::InfectedSevere] +
                                      last_value[(int)mio::osecir::InfectionState::InfectedCritical];
                if ((num_infected > switch_threshold * total_pop) || (num_infected < 1e-14)) {
                    return true;
                }
            }
            return false;
        };

        auto sim = HybridSim(std::move(sim_abm), std::move(sim_pbm), result_fct_abm, result_fct_pbm, true, params::t0,
                             params::dt_switch);
        timer.stop();
        init_time[run] = timer.get_elapsed_time();
        timer.reset();
        //Simulation
        timer.start();
        sim.advance(params::tmax, condition);
        timer.stop();
        sim_time[run]           = timer.get_elapsed_time();
        auto merged_result      = mio::merge_time_series(sim.get_result_model1(), sim.get_result_model2()).value();
        ensemble_result[run][0] = mio::interpolate_simulation_result(merged_result);
#pragma omp critical
        {
            if (sim.get_result_model2().get_num_time_points() > 1) {
                switch_runs++;
            }
        }
    }

    // Save ensemble results
    BOOST_OUTCOME_TRY(save_results(ensemble_result, result_dir, "hybrid"));
    mio::TimeSeries<double> init_time_ts(1);
    mio::TimeSeries<double> sim_time_ts(1);
    for (size_t i = 0; i < num_runs; i++) {
        init_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, init_time[i]));
        sim_time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, sim_time[i]));
    }
    BOOST_OUTCOME_TRY(init_time_ts.export_csv(result_dir + "hybrid_init_time.csv"));
    BOOST_OUTCOME_TRY(sim_time_ts.export_csv(result_dir + "hybrid_sim_time.csv"));
    std::cout << "Number of runs where switch occurred: " << switch_runs << " out of " << num_runs << std::endl;
    return mio::success();
}

int main()
{
    std::string base_dir                  = "/hpc_data/bick_ju/MemilioPaper/ScalingHybrid/";
    size_t num_runs                       = 100;
    std::vector<double> switch_thresholds = {0.02, 0.05, 0.1};
    std::vector<int> populations          = {100, 1000, 5000, 10000, 20000, 50000, 100000};
    double init_E                         = 0.01;

    mio::set_log_level(mio::LogLevel::err);

    for (auto pop : populations) {
        const std::string run_dir = base_dir + "pop_" + std::to_string(pop) + "_";
        auto abm_result           = simulate_dabm(run_dir, num_runs, pop, init_E);
        if (!abm_result) {
            printf("%s", "Error in dABM simulation\n");
            printf("%s\n", abm_result.error().formatted_message().c_str());
            return -1;
        }
        auto pbm_result = simulate_ode_secir(run_dir, num_runs, pop, init_E);
        if (!pbm_result) {
            printf("%s", "Error in OSECIR simulation\n");
            printf("%s\n", pbm_result.error().formatted_message().c_str());
            return -1;
        }
        for (auto threshold : switch_thresholds) {
            const std::string hybrid_dir = run_dir + std::to_string(static_cast<int>(threshold * 100)) + "_";
            auto hybrid_result           = simulate_hybrid(hybrid_dir, num_runs, threshold, pop, init_E);
            if (!hybrid_result) {
                printf("%s, %.02f", "Error in hybrid simulation with threshold\n", threshold);
                printf("%s\n", hybrid_result.error().formatted_message().c_str());
                return -1;
            }
        }
    }
    return 0;
}
