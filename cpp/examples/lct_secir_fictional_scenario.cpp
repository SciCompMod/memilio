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

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"

#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>

/** 
* @brief Perform a fictive simulation with realistic parameters and contacts, such that the reproduction number 
*   is approximately 1 at the beginning and rising or dropping at simulationtime 2.
*   
*   This scenario should enable a comparison of the qualitative behavior of different models.
*   
* @param[in] R0 Define R0 from simulationtime 2 on. Please use a number > 0.
* @param[in] num_subcompartments Number of subcompartments for each compartment where subcompartments make sense. 
*        Number is also used for some distributions of the IDE model.
*        Set num_subcompartments = 0 to use subcompartments with an expected sojourn time of approximately 1.
* @param[in] simulate_lct Defines if a simulation with an LCT model is done.
* @param[in] simulate_ide Defines if a simulation with an IDE model is done.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @param[in] tmax End time of the simulation.
* @param[in] print_result Specifies if the results should be printed.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate(ScalarType R0, int num_subcompartments = 3, bool simulate_lct = true,
                             bool simulate_ide = true, std::string save_dir = "", ScalarType tmax = 10,
                             bool print_result = false)
{

    ScalarType dt_flows         = 0.1;
    ScalarType total_population = 83155031.0;

    // Define parameters used for simulation and initialization.
    mio::lsecir::Parameters parameters_lct;
    parameters_lct.get<mio::lsecir::TimeExposed>()                      = 3.335;
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.31331;
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.94547;
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>()               = 11.634346;
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>()             = 17.476959;
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.0733271;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    if (R0 <= 1.) {
        // Perform simulation with dropping R0.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.7463));
        contact_matrix[0].add_damping(0., mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(R0, mio::SimulationTime(2.));
    }
    else {
        // Perform simulation with rising R0.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, R0 * 2.7463));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
    }

    contact_matrix[0].finalize();
    parameters_lct.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.3;
    parameters_lct.get<mio::lsecir::Seasonality>()                    = 0.;
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.206901;
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.0786429;
    parameters_lct.get<mio::lsecir::CriticalPerSevere>()              = 0.173176;
    parameters_lct.get<mio::lsecir::DeathsPerCritical>()              = 0.217177;

    // The initialization vector for the LCT model is calculated by defining transitions.
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
    init_transitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        SusceptibleToExposed_const * (1 - parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>());
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        (1 - parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>());
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        parameters_lct.get<mio::lsecir::CriticalPerSevere>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        (1 - parameters_lct.get<mio::lsecir::CriticalPerSevere>());
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        parameters_lct.get<mio::lsecir::DeathsPerCritical>();
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        (1 - parameters_lct.get<mio::lsecir::DeathsPerCritical>());
    init_transitions = init_transitions * dt_flows;

    // Add initial time point to time series.
    init.add_time_point(-350, init_transitions);
    // Add further time points until time 0 with constant values.
    while (init.get_last_time() < -1e-10) {
        init.add_time_point(init.get_last_time() + dt_flows, init_transitions);
    }

    // Set vector that specifies the number of subcompartments.
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    if (!(num_subcompartments == 0)) {
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = num_subcompartments;
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = num_subcompartments;
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = num_subcompartments;
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = num_subcompartments;
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = num_subcompartments;
    }
    else {
        // For approximately soujourn time of one day in each Subcompartment
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed] =
            round(parameters_lct.get<mio::lsecir::TimeExposed>());
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] =
            round(parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>());
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] =
            round(parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>());
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere] =
            round(parameters_lct.get<mio::lsecir::TimeInfectedSevere>());
        vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical] =
            round(parameters_lct.get<mio::lsecir::TimeInfectedCritical>());
    }
    mio::lsecir::InfectionState infectionState(vec_subcompartments);

    if (simulate_lct) {
        // Get initialization vector for LCT model with num_subcompartments subcompartments.
        mio::TimeSeries<ScalarType> init_copy(init);
        mio::lsecir::Initializer initializer(std::move(init_copy), infectionState, std::move(parameters_lct));
        Eigen::VectorXd init_compartments =
            initializer.compute_initializationvector(total_population, deaths, total_confirmed_cases);

        // Initialize model and perform simulation.
        mio::lsecir::Model model(std::move(init_compartments), infectionState, std::move(parameters_lct));
        mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(
            0, tmax, 0.5, model,
            std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(
                1e-10, 1e-5, 0, 0.1));
        // Calculate result without division in subcompartments.
        mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);

        if (print_result) {
            std::cout << "Result LCT model:" << std::endl;
            mio::lsecir::print_TimeSeries(populations, model.get_heading_CompartmentsBase());
        }
        if (!save_dir.empty()) {
            std::string R0string = std::to_string(R0);
            std::string filename = save_dir + "fictional_lct_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                                   std::to_string(num_subcompartments);
            if (tmax > 50) {
                filename = filename + "_long";
            }
            filename                               = filename + ".h5";
            mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
        }
    }

    //--------- IDE model ------------
    if (simulate_ide) {
        // Initialize model.
        mio::isecir::Model model_ide(std::move(init), total_population, deaths, total_confirmed_cases);

        // Set working parameters.
        // Set TransitionDistributions.
        mio::ConstantFunction initialfunc(0);
        mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
        std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistributioninit);

        mio::ExponentialDecay expdecayInfectedSevereToInfectedCritical(1. / 9.36);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]
            .set_state_age_function(expdecayInfectedSevereToInfectedCritical);

        mio::LognormSurvivalFunction lognInfectedSevereToRecovered(0.76, -0.45, 9.41);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_state_age_function(
            lognInfectedSevereToRecovered);

        mio::ExponentialDecay expdecayInfectedCriticalToDeath(1. / 14.88, 1);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_state_age_function(
            expdecayInfectedCriticalToDeath);
        expdecayInfectedCriticalToDeath.set_parameter(1 / 16.92);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_state_age_function(
            expdecayInfectedCriticalToDeath);

        mio::GammaSurvivalFunction erlangExposedToInfectedNoSymptoms(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed], 0,
            3.335 / vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
            erlangExposedToInfectedNoSymptoms);

        mio::GammaSurvivalFunction erlangInfectedNoSymptomsToInfectedSymptoms(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms], 0,
            1.865 / vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
            .set_state_age_function(erlangInfectedNoSymptomsToInfectedSymptoms);
        erlangInfectedNoSymptomsToInfectedSymptoms.set_scale(
            8.865 / vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
            erlangInfectedNoSymptomsToInfectedSymptoms);

        mio::GammaSurvivalFunction erlangInfectedSymptomsToInfectedSevere(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms], 0,
            6.30662 / vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
            .set_state_age_function(erlangInfectedSymptomsToInfectedSevere);
        erlangInfectedSymptomsToInfectedSevere.set_scale(
            7. / vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_state_age_function(
            erlangInfectedSymptomsToInfectedSevere);

        model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

        // Set other parameters.
        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            1 - parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            1 - parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            parameters_lct.get<mio::lsecir::CriticalPerSevere>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            1 - parameters_lct.get<mio::lsecir::CriticalPerSevere>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            parameters_lct.get<mio::lsecir::DeathsPerCritical>();
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            1 - parameters_lct.get<mio::lsecir::DeathsPerCritical>();

        model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

        model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        mio::ConstantFunction constfunc(parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>());
        mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
        model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
        StateAgeFunctionWrapperide.set_parameter(parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>());
        model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
        StateAgeFunctionWrapperide.set_parameter(parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>());
        model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

        model_ide.set_tol_for_support_max(1e-6);
        model_ide.check_constraints(dt_flows);

        // Simulate.
        mio::isecir::Simulation sim(model_ide, 0, dt_flows);
        sim.advance(tmax);

        if (print_result) {
            std::cout << "\nResult IDE model:" << std::endl;
            sim.print_compartments();
        }
        if (!save_dir.empty()) {
            std::string R0string     = std::to_string(R0);
            std::string filename_ide = save_dir + "fictional_ide_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                                       std::to_string(num_subcompartments);
            if (tmax > 50) {
                filename_ide = filename_ide + "_long";
            }
            filename_ide                           = filename_ide + ".h5";
            mio::IOResult<void> save_result_status = mio::save_result({sim.get_result()}, {0}, 1, filename_ide);
        }
    }
    return mio::success();
}

int main()
{
    // First case: drop R0 to 0.5.
    // Paths are valid if file is executed eg in memilio/build/bin
    std::string save_dir = "../../data/simulation_lct/dropR0/";
    auto result          = simulate(0.5, 0, true, true, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    result = simulate(0.5, 10, true, true, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    for (int i : {1, 3, 20}) {
        result = simulate(0.5, i, true, false, save_dir);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }

    // Second case: Rise R0 to 2 or 4.
    save_dir = "../../data/simulation_lct/riseR0short/";
    for (int j : {2, 4}) {
        for (int i : {0, 10}) {
            result = simulate(j, i, true, true, save_dir);
            if (!result) {
                printf("%s\n", result.error().formatted_message().c_str());
                return -1;
            }
        }
    }
    for (int j : {2, 4}) {
        for (int i : {1, 3, 20}) {
            result = simulate(j, i, true, false, save_dir);
            if (!result) {
                printf("%s\n", result.error().formatted_message().c_str());
                return -1;
            }
        }
    }
    // Third case: Rise R0 to 2 or 4 long term.
    save_dir = "../../data/simulation_lct/riseR0long/";

    for (int i : {0, 10}) {
        result = simulate(2, i, true, true, save_dir, 150);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }
    for (int i : {1, 3, 20}) {
        result = simulate(2, i, true, false, save_dir, 150);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }

    for (int i : {0, 10}) {
        result = simulate(4, i, true, true, save_dir, 75);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }
    for (int i : {1, 3, 20}) {
        result = simulate(4, i, true, false, save_dir, 75);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }
    return 0;
}