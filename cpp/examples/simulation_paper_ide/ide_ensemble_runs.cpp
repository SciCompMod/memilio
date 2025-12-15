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
#include "models/ide_secir/model.h"
#include "models/ide_secir/infection_state.h"
#include "models/ide_secir/simulation.h"
#include "models/ide_secir/parameters_io.h"

#include "memilio/compartments/parameter_studies.h"
#include "memilio/utils/base_dir.h"
#include "memilio/io/cli.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"
#include "memilio/utils/miompi.h"

#include <omp.h>
#include <mpi.h>
#include "boost/filesystem.hpp"
#include <string>

ScalarType uncertain(ScalarType v)
{
    const double var = .1;
    return mio::ParameterDistributionUniform(v * (1 - var), v * (1 + var)).get_sample(mio::thread_local_rng());
}

namespace params
{
size_t num_agegroups = 6;
size_t total_population = 80 * 1e6;
int num_processes = 1;

// Epidemiological parameters

// Define (age-resolved) parameters.
mio::Date start_date(2020, 12, 24);
const ScalarType t0                       = 0;
const ScalarType tmax                     = 30;
const ScalarType dt                       = 0.01;

}

using Model = mio::isecir::Model;

mio::IOResult<Model> initialize_isecir(std::string data_dir) 
{
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(params::num_agegroups));
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_init = mio::CustomIndexArray<ScalarType, mio::AgeGroup>(
        mio::AgeGroup(params::num_agegroups),
        0.); // The number of deaths will be overwritten if reported data is used for initialization.
    
    Model model(mio::TimeSeries<ScalarType>((size_t)mio::isecir::InfectionTransition::Count),
                             total_population_init, deaths_init, params::num_agegroups);

    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(mio::path_join(data_dir, "Germany/pydata/cases_all_age_ma7.json")));
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> scale_confirmed_cases =
            mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(params::num_agegroups), 1.);
    BOOST_OUTCOME_TRY(mio::isecir::set_initial_flows<mio::ConfirmedCasesDataEntry>(
            model, params::dt, rki_data, params::start_date, scale_confirmed_cases));

    return mio::success(model);
}

Model draw_sample(const Model& model)
{
    auto copy = model;

    // TransitionDistributions
    mio::SmootherCosine<ScalarType> smoothcos1(uncertain(2.0));
    mio::StateAgeFunctionWrapper<ScalarType> delaydistribution1(smoothcos1);
    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib1((size_t)mio::isecir::InfectionTransition::Count, delaydistribution1);
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(params::num_agegroups); ++group) {
        copy.parameters.get<mio::isecir::TransitionDistributions>()[group] = vec_delaydistrib1;
    }

    // Furhter epidemiological parameters
    auto draw_const_func = [&](ScalarType value){
        mio::ConstantFunction<ScalarType> constfunc(uncertain(value));
        return mio::StateAgeFunctionWrapper<ScalarType>(constfunc);
    };
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(params::num_agegroups); ++group) {
        copy.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = draw_const_func(1.0);
        copy.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = draw_const_func(1.0);
        copy.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = draw_const_func(1.0);
    }

    return copy;
}

mio::IOResult<void> simulate(std::string save_dir, std::string data_dir, size_t num_ensemble_runs)
{
    mio::set_log_level(mio::LogLevel::off);

    using namespace params;
    if (mio::mpi::is_root()) {
        std::cout << "Realistic scenario." << std::endl;
    }

    Model model = initialize_isecir(data_dir).value();
    // Set integrator of fifth order with fixed step size and perform simulation.
    // auto integrator =
    //     std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // // Choose dt_min = dt_max to get a fixed step size.
    // integrator->set_dt_min(dt);
    // integrator->set_dt_max(dt);

    mio::ParameterStudy parameter_study(model, params::t0, params::tmax, params::dt, num_ensemble_runs);
    ScalarType total_time = 0;
    if (mio::mpi::is_root()) {
        total_time -= omp_get_wtime();
    }
    parameter_study.run(
        [](auto&& params_model, ScalarType, ScalarType dt, size_t) {
            auto copy = params_model;
            return mio::isecir::Simulation(draw_sample(copy), dt);
        }, 
        [&](auto results_model, auto&& run_idx) {
            auto interpolated_results = mio::interpolate_simulation_result(results_model.get_result());
            mio::unused(interpolated_results, run_idx);
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
    mio::set_log_level(mio::LogLevel::warn);
    mio::unused(save_dir);
    return mio::success();
}

int main(int argc, char** argv)
{
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                          .add<"ResultDirectory">(mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper_ide/results_ensemble"))
                          .add<"DataDirectory">(mio::path_join(mio::base_dir(), "data"))
                          .add<"NumberEnsembleRuns">(100, {.alias = "nRun"})
                          .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"ResultDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }

    boost::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    boost::filesystem::create_directories(res_dir);

    mio::mpi::init();
    int size;
    MPI_Comm_size(mio::mpi::get_world(), &size);
    params::num_processes = size;

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