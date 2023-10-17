/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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
#include "lct_secir/initialization.h"
#include "lct_secir/parameters.h"
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>

/** 
* @brief Perform simulation for an LCT model with realistic parameters with different initialization methods to compare those.
*
* @param[in] num_subcompartments Number of subcompartments for each compartment where subcompartments make sense.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @param[in] tmax End time of the simulation.
* @param[in] print_result Specifies if the results should be printed.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate(int num_subcompartments = 3, std::string save_dir = "", ScalarType tmax = 20,
                             bool print_result = false)
{
    ScalarType dt_flows         = 0.1;
    ScalarType total_population = 83155031.0;

    // Define parameters used for simulation and initialization.
    // Parameters are calculated via examples/compute_parameters.cpp.
    mio::lsecir::Parameters parameters;
    parameters.get<mio::lsecir::TimeExposed>()                      = 3.335;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.31331;
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.94547;
    parameters.get<mio::lsecir::TimeInfectedSevere>()               = 11.634346;
    parameters.get<mio::lsecir::TimeInfectedCritical>()             = 17.476959;
    parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.0733271;

    mio::ContactMatrixGroup contact_matrix         = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                              = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.7463));
    parameters.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.3;
    parameters.get<mio::lsecir::Seasonality>()                    = 0;
    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.206901;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.0786429;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.173176;
    parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.217177;

    // --- First initialization method: The initialization vector for the LCT model is calculated by defining transitions.
    // Create TimeSeries with num_transitions elements.
    int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization of transitions.
    /* For this example, the intention is to create nearly constant values for SusceptiblesToExposed flow 
    at the beginning of the simulation. Therefore we initalize the flows accordingly constant for 
    SusceptiblesToExposed and derive matching values for the other flows.*/
    // 7-Tage-Inzidenz at 15.10.2020 was 34.1, see https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Okt_2020/2020-10-15-de.pdf?__blob=publicationFile.
    ScalarType SusceptibleToExposed_const = (34.1 / 7) * total_population / 100000;
    ScalarType total_confirmed_cases      = 341223;
    ScalarType deaths                     = 9710;
    Eigen::VectorXd init_transitions(num_transitions);
    init_transitions[(int)mio::lsecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
    init_transitions[(int)mio::lsecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        SusceptibleToExposed_const * (1 - parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>());
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>();
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        parameters.get<mio::lsecir::SeverePerInfectedSymptoms>();
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToRecovered] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        (1 - parameters.get<mio::lsecir::SeverePerInfectedSymptoms>());
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        parameters.get<mio::lsecir::CriticalPerSevere>();
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToRecovered] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        (1 - parameters.get<mio::lsecir::CriticalPerSevere>());
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedCriticalToDead] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        parameters.get<mio::lsecir::DeathsPerCritical>();
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedCriticalToRecovered] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        (1 - parameters.get<mio::lsecir::DeathsPerCritical>());
    init_transitions = init_transitions * dt_flows;

    // Add initial time point to time series.
    init.add_time_point(-200, init_transitions);
    // Add further time points until time 0.
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt_flows, init_transitions);
    }

    // Set vector that specifies the number of subcompartments.
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = num_subcompartments;
    mio::lsecir::InfectionState infectionState(vec_subcompartments);

    mio::TimeSeries<ScalarType> init_copy(init);

    // Get initialization vector for LCT model with num_subcompartments subcompartments.
    mio::lsecir::Initializer initializer(std::move(init_copy), infectionState, std::move(parameters));
    auto init_compartments = initializer.compute_initializationvector(total_population, deaths, total_confirmed_cases);

    // Initialize model and perform simulation.
    mio::lsecir::Model model(std::move(init_compartments), infectionState, std::move(parameters));
    mio::TimeSeries<ScalarType> result_transitions = mio::lsecir::simulate(
        0, tmax, 0.5, model,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(1e-10, 1e-5, 0,
                                                                                                         0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations_transitions = model.calculate_populations(result_transitions);

    if (print_result) {
        mio::lsecir::print_TimeSeries(populations_transitions, model.get_heading_CompartmentsBase());
    }
    if (!save_dir.empty()) {
        std::string filename = save_dir + "lct_init_transitions_" + std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename = filename + "_long";
        }
        filename                               = filename + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({populations_transitions}, {0}, 1, filename);
    }

    // --- Second initialization method: The initialization vector for the LCT model is calculated by taking expected sojourn times.
    // Because of the constant initialization, the constant value of SusceptibleToExposed is scaled by expected sojourn time.
    Eigen::VectorXd init_vec((int)mio::lsecir::InfectionStateBase::Count);
    init_vec[(int)mio::lsecir::InfectionStateBase::Exposed] =
        SusceptibleToExposed_const * parameters.get<mio::lsecir::TimeExposed>();
    init_vec[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] =
        SusceptibleToExposed_const * parameters.get<mio::lsecir::TimeInfectedNoSymptoms>();
    init_vec[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] =
        SusceptibleToExposed_const * (1 - parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()) *
        parameters.get<mio::lsecir::TimeInfectedSymptoms>();
    init_vec[(int)mio::lsecir::InfectionStateBase::InfectedSevere] =
        SusceptibleToExposed_const * (1 - parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()) *
        parameters.get<mio::lsecir::SeverePerInfectedSymptoms>() * parameters.get<mio::lsecir::TimeInfectedSevere>();
    init_vec[(int)mio::lsecir::InfectionStateBase::InfectedCritical] =
        SusceptibleToExposed_const * (1 - parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()) *
        parameters.get<mio::lsecir::SeverePerInfectedSymptoms>() * parameters.get<mio::lsecir::CriticalPerSevere>() *
        parameters.get<mio::lsecir::TimeInfectedCritical>();
    init_vec[(int)mio::lsecir::InfectionStateBase::Recovered] =
        total_confirmed_cases - deaths - init_vec[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] -
        init_vec[(int)mio::lsecir::InfectionStateBase::InfectedSevere] -
        init_vec[(int)mio::lsecir::InfectionStateBase::InfectedCritical];
    init_vec[(int)mio::lsecir::InfectionStateBase::Dead] = deaths;
    init_vec[(int)mio::lsecir::InfectionStateBase::Susceptible] =
        total_population - init_vec[(int)mio::lsecir::InfectionStateBase::Exposed] -
        init_vec[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] -
        init_vec[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] -
        init_vec[(int)mio::lsecir::InfectionStateBase::InfectedSevere] -
        init_vec[(int)mio::lsecir::InfectionStateBase::InfectedCritical] -
        init_vec[(int)mio::lsecir::InfectionStateBase::Recovered] -
        init_vec[(int)mio::lsecir::InfectionStateBase::Dead];

    // Constant initialization leads to equally distributed values accross substates.
    Eigen::VectorXd init_mean(infectionState.get_count());
    for (int i = 0; i < (int)mio::lsecir::InfectionStateBase::Count; i++) {
        for (int j = infectionState.get_firstindex(i);
             j < infectionState.get_firstindex(i) + infectionState.get_number(i); j++) {
            init_mean[j] = init_vec[i] / infectionState.get_number(i);
        }
    }

    // Initialize model and perform simulation.
    mio::lsecir::Model model2(std::move(init_mean), infectionState, std::move(parameters));
    mio::TimeSeries<ScalarType> result2 = mio::lsecir::simulate(
        0, tmax, 0.5, model2,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(1e-10, 1e-5, 0,
                                                                                                         0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations2 = model2.calculate_populations(result2);

    if (print_result) {
        mio::lsecir::print_TimeSeries(populations2, model2.get_heading_CompartmentsBase());
    }
    if (!save_dir.empty()) {
        std::string filename_mean = save_dir + "lct_init_mean_" + std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename_mean = filename_mean + "_long";
        }
        filename_mean                           = filename_mean + ".h5";
        mio::IOResult<void> save_result_status2 = mio::save_result({populations2}, {0}, 1, filename_mean);
    }

    // --- Third initialization method: Initial values are distributed only in the first substates.
    Eigen::VectorXd init_first = Eigen::VectorXd::Zero(infectionState.get_count());
    for (int i = 0; i < (int)mio::lsecir::InfectionStateBase::Count; i++) {
        init_first[infectionState.get_firstindex(i)] = init_vec[i];
    }

    // Initialize model and perform simulation.
    mio::lsecir::Model model3(std::move(init_first), infectionState, std::move(parameters));
    mio::TimeSeries<ScalarType> result3 = mio::lsecir::simulate(
        0, tmax, 0.5, model3,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(1e-10, 1e-5, 0,
                                                                                                         0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations3 = model3.calculate_populations(result3);

    if (print_result) {
        mio::lsecir::print_TimeSeries(populations3, model3.get_heading_CompartmentsBase());
    }
    if (!save_dir.empty()) {
        std::string filename_first = save_dir + "lct_init_first_" + std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename_first = filename_first + "_long";
        }
        filename_first                          = filename_first + ".h5";
        mio::IOResult<void> save_result_status3 = mio::save_result({populations3}, {0}, 1, filename_first);
    }
    return mio::success();
}

int main()
{ // Path is valid if file is executed eg in memilio/build/bin
    std::string save_dir = "../../data/simulation_lct/init/";
    auto result          = simulate(20, save_dir, 100);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    for (int i : {3, 10, 20}) {
        result = simulate(i, save_dir, 20);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }
    return 0;
}