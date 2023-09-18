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
#include "memilio/utils/time_series.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <math.h>
#include <iostream>

int main()
{
    bool save_result   = true;
    bool print_result  = false;
    bool simulate_ide  = true;
    bool simulate_lct  = false;
    using Vec          = mio::TimeSeries<ScalarType>::Vector;
    using ParameterSet = mio::lsecir::Parameters;

    ScalarType dt_flows         = 0.1;
    ScalarType total_population = 83155031.0;
    ScalarType tmax             = 20;

    // Define ParameterSet used for simulation and initialization.
    ParameterSet parameters_lct;
    parameters_lct.get<mio::lsecir::TimeExposed>()                      = 3.335;
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.31331;
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.94547;
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>()               = 11.634346;
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>()             = 17.476959;
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.0733271;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    // Perform simulation with dropping R0.
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.7463));
    contact_matrix[0].add_damping(0., mio::SimulationTime(4.9));
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(5.));
    // Perform simulation with rising R0.
    /*contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2 * 2.7463));
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(-200.));
    contact_matrix[0].add_damping(0.5, mio::SimulationTime(4.9));
    contact_matrix[0].add_damping(0., mio::SimulationTime(5.));*/

    contact_matrix[0].finalize();
    parameters_lct.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.3;
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
    Vec init_transitions(num_transitions);
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
    /*int num_subcompartments = 20;
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = num_subcompartments;
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = num_subcompartments;
    mio::lsecir::InfectionState infectionStates(vec_subcompartments);*/

    // For approximately soujourn time of one day in each Subcompartment
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
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
    mio::lsecir::InfectionState infectionStates(vec_subcompartments);
    if (simulate_lct) {
        // Get initialization vector for LCT model with num_subcompartments subcompartments.
        mio::TimeSeries<ScalarType> init_copy(init);
        mio::lsecir::Initializer initializer(std::move(init_copy), infectionStates, std::move(parameters_lct));
        auto init_compartments =
            initializer.compute_initializationvector(total_population, deaths, total_confirmed_cases);

        // Initialize model and perform simulation.
        mio::lsecir::Model model(std::move(init_compartments), infectionStates, std::move(parameters_lct));
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
        if (save_result) {
            //auto save_result_status_subcompartments =
            //   mio::save_result({result}, {0}, 1, "result_lct_subcompartments_fictional_1.h5");
            auto save_result_status = mio::save_result({populations}, {0}, 1, "result_lct_fictional_var.h5");
        }
    }

    //--------- IDE model ------------
    // Initialize model.
    if (simulate_ide) {
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
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed] / 3.335,
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
            erlangExposedToInfectedNoSymptoms);

        mio::GammaSurvivalFunction erlangInfectedNoSymptomsToInfectedSymptoms(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] / 1.865,
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
            .set_state_age_function(erlangInfectedNoSymptomsToInfectedSymptoms);
        erlangInfectedNoSymptomsToInfectedSymptoms.set_rate(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] / 8.865);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
            erlangInfectedNoSymptomsToInfectedSymptoms);

        mio::GammaSurvivalFunction erlangInfectedSymptomsToInfectedSevere(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] / 7.64507,
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
            .set_state_age_function(erlangInfectedSymptomsToInfectedSevere);
        erlangInfectedSymptomsToInfectedSevere.set_rate(
            vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] / 7.);
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

        // Carry out simulation.
        mio::isecir::Simulation sim(model_ide, 0, dt_flows);
        sim.advance(tmax);

        if (print_result) {
            std::cout << "\nResult IDE model:" << std::endl;
            sim.print_compartments();
        }
        if (save_result) {
            auto save_result_status = mio::save_result({sim.get_result()}, {0}, 1, "result_ide_fictional.h5");
        }
    }
}