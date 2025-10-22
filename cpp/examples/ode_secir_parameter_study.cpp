/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/compartments/simulation.h"
#include "memilio/config.h"
#include "ode_secir/model.h"
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
mio::IOResult<void> write_single_run_result(const size_t run, const mio::osecir::Simulation<ScalarType>& sim)
{
    std::string abs_path;
    BOOST_OUTCOME_TRY(auto&& created, mio::create_directory("results", abs_path));

    if (run == 0) {
        std::cout << "Results are stored in " << abs_path << '\n';
        if (!created) {
            std::cout << "Directory already exists, files from previous runs will be overwritten." << '\n';
        }
    }

    //write sampled parameters for this run
    auto node_filename = mio::path_join(abs_path, "Parameters_run" + std::to_string(run) + ".json");
    BOOST_OUTCOME_TRY(mio::write_json(node_filename, sim.get_result()));

    //write results for this run
    std::vector<mio::TimeSeries<ScalarType>> all_results;
    std::vector<int> ids;

    BOOST_OUTCOME_TRY(mio::save_result({sim.get_result()}, {0}, (int)sim.get_model().parameters.get_num_groups().get(),
                                       mio::path_join(abs_path, ("Results_run" + std::to_string(run) + ".h5"))));

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 50;
    ScalarType dt   = 0.1;

    ScalarType cont_freq = 10; // see Polymod study

    ScalarType num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20,
               num_icu_t0 = 10, num_rec_t0 = 10, num_dead_t0 = 0;

    mio::osecir::Model<ScalarType> model(1);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    ScalarType fact          = 1.0 / (ScalarType)(size_t)num_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity<ScalarType>>(std::numeric_limits<ScalarType>::max());
    params.set<mio::osecir::StartDay<ScalarType>>(0);
    params.set<mio::osecir::Seasonality<ScalarType>>(0);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::osecir::TimeExposed<ScalarType>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 2.;
        params.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 6.;
        params.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[i]     = 12;
        params.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[i]   = 8;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<ScalarType>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<ScalarType>>()[i]                = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup<ScalarType>& contact_matrix = params.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                                   = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal<ScalarType>(model, t0, tmax, 0.2);

    auto write_parameters_status = mio::write_json("parameters.json", model);
    if (!write_parameters_status) {
        std::cout << "Error writing parameters: " << write_parameters_status.error().formatted_message();
        return -1;
    }

    //create study
    auto num_runs = size_t(1);
    // mio::ParameterStudy2<mio::osecir::Simulation<ScalarType>, mio::osecir::Model<ScalarType>, ScalarType>
    //     parameter_study(model, t0, tmax, dt, num_runs);

    auto parameter_study =
        mio::make_parameter_study<mio::osecir::Simulation<ScalarType>>(model, t0, tmax, dt, num_runs);

    //run study
    auto sample_graph = [](const auto& model_, ScalarType t0_, ScalarType dt_, size_t) {
        mio::osecir::Model<ScalarType> copy = model_;
        mio::osecir::draw_sample(copy);
        return mio::osecir::Simulation<ScalarType>(std::move(copy), t0_, dt_);
    };
    auto handle_result = [](auto&& sim, auto&& run) {
        auto write_result_status = write_single_run_result(run, sim);
        if (!write_result_status) {
            std::cout << "Error writing result: " << write_result_status.error().formatted_message();
        }
        return 0; //Result handler must return something, but only meaningful when using MPI.
    };
    parameter_study.run(sample_graph, handle_result);

    return 0;
}
