/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"

#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"

#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "ode_secir/infection_state.h"
#include <string>
#include <map>
#include <iostream>

using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

// necessary because num_subcompartments is used as a template argument and has ti be
constexpr int num_subcompartments = 20;

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"t0", 0.},
                                                          {"dt_flows", 0.1},
                                                          {"total_population", 83155031.},
                                                          {"total_confirmed_cases", 341223.},
                                                          {"deaths", 9710.},
                                                          {"TimeExposed", 3.335},
                                                          {"TimeInfectedNoSymptoms", 3.31331},
                                                          {"TimeInfectedSymptoms", 6.94547},
                                                          {"TimeInfectedSevere", 11.634346},
                                                          {"TimeInfectedCritical", 17.476959},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"RecoveredPerInfectedNoSymptoms", 0.206901},
                                                          {"SeverePerInfectedSymptoms", 0.0786429},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.217177}};

mio::UncertainContactMatrix<ScalarType> get_contact_matrix(ScalarType R0)
{
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

    return mio::UncertainContactMatrix(contact_matrix);
}

mio::TimeSeries<ScalarType> get_initial_flows()
{
    // The initialization vector for the LCT model is calculated by defining transitions.
    // Create TimeSeries with num_transitions elements.
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization of transitions.
    /* For this example, the intention is to create nearly constant values for SusceptiblesToExposed flow
    at the beginning of the simulation. Therefore we initalize the flows accordingly constant for
    SusceptiblesToExposed and derive matching values for the other flows.*/
    // 7-Tage-Inzidenz at 15.10.2020 was 34.1, see https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Okt_2020/2020-10-15-de.pdf?__blob=publicationFile.
    ScalarType SusceptibleToExposed_const = (34.1 / 7.) * simulation_parameter["total_population"] / 100000.;
    Eigen::VectorXd init_transitions(num_transitions);
    init_transitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        SusceptibleToExposed_const * (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        simulation_parameter["SeverePerInfectedSymptoms"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        (1 - simulation_parameter["SeverePerInfectedSymptoms"]);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        simulation_parameter["CriticalPerSevere"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        (1 - simulation_parameter["CriticalPerSevere"]);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        simulation_parameter["DeathsPerCritical"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        (1 - simulation_parameter["DeathsPerCritical"]);
    init_transitions = init_transitions * simulation_parameter["dt_flows"];

    // Add initial time point to time series.
    init.add_time_point(-350, init_transitions);
    // Add further time points until time 0 with constant values.
    while (init.get_last_time() < simulation_parameter["t0"] - 1e-10) {
        init.add_time_point(init.get_last_time() + simulation_parameter["dt_flows"], init_transitions);
    }
    return init;
}

mio::TimeSeries<ScalarType> simulate_ide_model(ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    // Initialize model.
    mio::isecir::Model model_ide(std::move(get_initial_flows()), simulation_parameter["total_population"],
                                 simulation_parameter["deaths"], simulation_parameter["total_confirmed_cases"]);

    // Set working parameters.
    // Set TransitionDistributions.
    mio::ConstantFunction initialfunc(0);
    mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistributioninit);

    mio::ExponentialSurvivalFunction expInfectedSevereToInfectedCritical(1. / 9.36);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical].set_state_age_function(
        expInfectedSevereToInfectedCritical);

    mio::LognormSurvivalFunction lognInfectedSevereToRecovered(0.76, -0.45, 9.41);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_state_age_function(
        lognInfectedSevereToRecovered);

    mio::ExponentialSurvivalFunction expInfectedCriticalToDeath(1. / 14.88, 1);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_state_age_function(
        expInfectedCriticalToDeath);
    expInfectedCriticalToDeath.set_distribution_parameter(1 / 16.92);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_state_age_function(
        expInfectedCriticalToDeath);

    mio::GammaSurvivalFunction erlangExposedToInfectedNoSymptoms(num_subcompartments, 0, 3.335 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
        erlangExposedToInfectedNoSymptoms);

    mio::GammaSurvivalFunction erlangInfectedNoSymptomsToInfectedSymptoms(num_subcompartments, 0,
                                                                          1.865 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_state_age_function(erlangInfectedNoSymptomsToInfectedSymptoms);
    erlangInfectedNoSymptomsToInfectedSymptoms.set_scale(8.865 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
        erlangInfectedNoSymptomsToInfectedSymptoms);

    mio::GammaSurvivalFunction erlangInfectedSymptomsToInfectedSevere(num_subcompartments, 0,
                                                                      6.30662 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
        erlangInfectedSymptomsToInfectedSevere);
    erlangInfectedSymptomsToInfectedSevere.set_scale(7. / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_state_age_function(
        erlangInfectedSymptomsToInfectedSevere);

    model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set other parameters.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
        1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
        1 - simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
        simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
        1 - simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
        simulation_parameter["DeathsPerCritical"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
        1 - simulation_parameter["DeathsPerCritical"];

    model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    model_ide.parameters.get<mio::isecir::ContactPatterns>() = get_contact_matrix(R0);

    mio::ConstantFunction constfunc(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

    model_ide.set_tol_for_support_max(1e-6);
    model_ide.check_constraints(simulation_parameter["dt_flows"]);

    // Simulate.
    mio::isecir::Simulation sim(model_ide, simulation_parameter["dt_flows"]);
    sim.advance(tmax);

    if (!save_dir.empty()) {

        std::string R0string     = std::to_string(R0);
        std::string filename_ide = save_dir + "fictional_ide_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                                   std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename_ide = filename_ide + "_long";
        }
        std::string filename_ide_flows = filename_ide + "_flows.h5";
        mio::IOResult<void> save_result_status_f =
            mio::save_result({sim.get_transitions()}, {0}, 1, filename_ide_flows);
        std::string filename_ide_compartments = filename_ide + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({sim.get_result()}, {0}, 1, filename_ide_compartments);
    }

    // Return vector with initial compartments.
    return sim.get_result();
}

void simulate_ode_model(Vector init_compartments, ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    // auto init_compartments = init_compartments2.get_value(0);
    // Use FlowModel to make results directly comparable to IDE model.
    mio::osecir::Model model_ode(1);

    // Set initial values for compartments.
    // Use mio::isecir::InfectionState when accessing init_compartments since this is computed using the IDE model.
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] =
        init_compartments[int(mio::isecir::InfectionState::Exposed)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedNoSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedSevere)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedCritical)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}] =
        init_compartments[int(mio::isecir::InfectionState::Recovered)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}] =
        init_compartments[int(mio::isecir::InfectionState::Dead)];
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    simulation_parameter["total_population"]);

    // Set working parameters.
    model_ode.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeExposed"];
    model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedNoSymptoms"];
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedSymptoms"];
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedSevere"];
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedCritical"];

    // Set probabilities that determine proportion between compartments.
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["CriticalPerSevere"];
    model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["DeathsPerCritical"];

    // Further model parameters.
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TransmissionProbabilityOnContact"];
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic.
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0.
    model_ode.parameters.get<mio::osecir::ICUCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix.
    model_ode.parameters.set<mio::osecir::Seasonality<ScalarType>>(simulation_parameter["Seasonality"]);

    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() = get_contact_matrix(R0);

    model_ode.check_constraints();

    // auto integrator =
    //     std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // integrator->set_dt_min(simulation_parameter["dt_flows"]);
    // integrator->set_dt_max(simulation_parameter["dt_flows"]);

    std::vector<mio::TimeSeries<ScalarType>> results_ode =
        simulate_flows(simulation_parameter["t0"], tmax, simulation_parameter["dt_flows"], model_ode);

    // flows_ode.print_table();

    if (!save_dir.empty()) {
        std::string R0string     = std::to_string(R0);
        std::string filename_ode = save_dir + "fictional_ode_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                                   std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename_ode = filename_ode + "_long";
        }
        std::string filename_ode_flows           = filename_ode + "_flows.h5";
        mio::IOResult<void> save_result_status_f = mio::save_result({results_ode[1]}, {0}, 1, filename_ode_flows);
        std::string filename_ode_compartments    = filename_ode + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({results_ode[0]}, {0}, 1, filename_ode_compartments);
    }
}

int main()
{
    // Paths are valid if file is executed eg in memilio/build/bin.
    std::string save_dir = "../../results/";
    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directory(dir);

    // Options used: For R0=2 epidemic peak use tmax=150,
    // for R0=4 epidemic peak use tmax = 75.
    // For short things: 10 days and R0=0.5 or 2
    ScalarType R0   = 0.5;
    ScalarType tmax = 10;

    mio::TimeSeries<ScalarType> result = simulate_ide_model(R0, tmax, save_dir);
    // if (!result) {
    //     printf("%s\n", result.error().formatted_message().c_str());
    //     return -1;
    // }

    Vector compartments = result.get_value(0);
    ;

    simulate_ode_model(compartments, R0, tmax, save_dir);

    R0   = 2.;
    tmax = 10;

    result = simulate_ide_model(R0, tmax, save_dir);
    // if (!result) {
    //     printf("%s\n", result.error().formatted_message().c_str());
    //     return -1;
    // }
    compartments = result.get_value(0);
    simulate_ode_model(compartments, R0, tmax, save_dir);

    return 0;
}