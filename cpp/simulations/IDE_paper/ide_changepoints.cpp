/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"

#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "ode_secir/infection_state.h"
#include <string>
#include <map>

using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

// Used parameters
std::map<std::string, ScalarType> simulation_parameter = {
    {"t0", 0.},
    {"dt", 0.01},
    {"total_population", 83155031.},
    {"total_confirmed_cases", 341223.},
    {"deaths", 0.},
    {"TimeExposed", 4.5},
    {"TimeInfectedNoSymptoms", 2.527617},
    {"TimeInfectedSymptoms", 7.889900},
    {"TimeInfectedSevere", 15.225278},
    {"TimeInfectedCritical", 15.230258},
    {"TransmissionProbabilityOnContact", 0.0733271},
    {"RelativeTransmissionNoSymptoms", 1},
    {"RiskOfInfectionFromSymptomatic", 0.3},
    {"Seasonality", 0.},
    {"InfectedSymptomsPerInfectedNoSymptoms", 0.793099},
    {"SeverePerInfectedSymptoms", 0.078643},
    {"CriticalPerSevere", 0.173176},
    {"DeathsPerCritical", 0.387803},
    {"cont_freq", 3.114219}}; // Computed so that we obtain constant new infections at beginning of simulation.

/**
* @brief Function to scale the contact matrix according to factor contact_scaling after two days.  
*
* @param[in] contact_scaling Factor that is applied to contact matrix after two days. 
* @returns Scaled contact matrix.
*/
mio::UncertainContactMatrix<ScalarType> scale_contact_matrix(ScalarType contact_scaling)
{
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    if (contact_scaling <= 1.) {
        // Perform simulation with a decrease in contacts.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, simulation_parameter["cont_freq"]));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
        contact_matrix[0].add_damping(contact_scaling, mio::SimulationTime(2.1));
    }
    else {
        // Perform simulation with an increase in contacts.
        contact_matrix[0] =
            mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, contact_scaling * simulation_parameter["cont_freq"]));
        contact_matrix[0].add_damping(1 - 1. / contact_scaling, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / contact_scaling, mio::SimulationTime(2.));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.1));
    }

    return mio::UncertainContactMatrix(contact_matrix);
}

/**
* @brief Function to compute the initial flows needed for the IDE model where we assume that we have a constant number 
* of new transmissions.  
*
* @returns TimeSeries containing intitial flows. 
*/
mio::TimeSeries<ScalarType> get_initial_flows()
{
    // The initialization vector for the IDE model is calculated by defining transitions.
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
        SusceptibleToExposed_const * simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * (1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"]);
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
    init_transitions = init_transitions * simulation_parameter["dt"];

    // Add initial time point to time series.
    init.add_time_point(-350, init_transitions);
    // Add further time points until time 0 with constant values.
    while (init.get_last_time() < simulation_parameter["t0"] - 1e-3) {
        init.add_time_point(init.get_last_time() + simulation_parameter["dt"], init_transitions);
    }
    return init;
}

/**
* @brief Function that simulates from time 0 until tmax using an IDE model where we apply a contact scaling after
* two days.   
*
* @param[in] contact_scaling Factor that is applied to contact matrix after two days. 
* @param[in] tmax Time up to which we simulate. 
* @param[in] save_dir Directory where simulation results will be stored. 
* @returns Any io errors that happen.
*/
mio::IOResult<mio::TimeSeries<ScalarType>> simulate_ide_model(ScalarType contact_scaling, ScalarType tmax,
                                                              std::string save_dir = "")
{
    // Initialize model.
    size_t num_agegroups = 1;
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups),
                                                         simulation_parameter["total_population"]);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), simulation_parameter["deaths"]);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups),
                                                         simulation_parameter["total_confirmed_cases"]);

    mio::isecir::Model model_ide(get_initial_flows(), total_population, deaths, num_agegroups, total_confirmed_cases);

    // Set working parameters.
    // Set TransitionDistributions.
    mio::ConstantFunction initialfunc(0);
    mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistributioninit);
    // ExposedToInfectedNoSymptoms
    mio::LognormSurvivalFunction survivalExposedToInfectedNoSymptoms(0.32459285, 0, 4.26907484);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
        survivalExposedToInfectedNoSymptoms);
    // InfectedNoSymptomsToInfectedSymptoms
    mio::LognormSurvivalFunction survivalInfectedNoSymptomsToInfectedSymptoms(0.71587510, 0, 0.85135303);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_state_age_function(survivalInfectedNoSymptomsToInfectedSymptoms);
    // InfectedNoSymptomsToRecovered
    mio::LognormSurvivalFunction survivalInfectedNoSymptomsToRecovered(0.24622068, 0, 7.7611400);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
        survivalInfectedNoSymptomsToRecovered);
    // InfectedSymptomsToInfectedSevere
    mio::LognormSurvivalFunction survivalInfectedSymptomsToInfectedSevere(0.66258947, 0, 5.29920733);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
        survivalInfectedSymptomsToInfectedSevere);
    // InfectedSymptomsToRecovered
    mio::LognormSurvivalFunction survivalInfectedSymptomsToRecovered(0.24622068, 0, 7.76114000);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_state_age_function(
        survivalInfectedSymptomsToRecovered);
    // InfectedSevereToInfectedCritical
    mio::LognormSurvivalFunction survivalInfectedSevereToInfectedCritical(1.01076765, 0, 0.90000000);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical].set_state_age_function(
        survivalInfectedSevereToInfectedCritical);
    // InfectedSevereToRecovered
    mio::LognormSurvivalFunction survivalInfectedSevereToRecovered(0.33816427, 0, 17.09411753);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_state_age_function(
        survivalInfectedSevereToRecovered);
    // InfectedCriticalToDead
    mio::LognormSurvivalFunction survivalInfectedCriticalToDead(0.42819924, 0, 9.76267505);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_state_age_function(
        survivalInfectedCriticalToDead);
    // InfectedCriticalToRecovered
    mio::LognormSurvivalFunction survivalInfectedCriticalToRecovered(0.33816427, 0, 17.09411753);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_state_age_function(
        survivalInfectedCriticalToRecovered);

    model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set other parameters.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
        simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
        1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
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

    model_ide.parameters.get<mio::isecir::ContactPatterns>() = scale_contact_matrix(contact_scaling);

    mio::ConstantFunction constfunc(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

    model_ide.set_tol_for_support_max(1e-6);
    model_ide.check_constraints(simulation_parameter["dt"]);

    // Simulate.
    mio::isecir::Simulation sim(model_ide, simulation_parameter["dt"]);
    sim.advance(tmax);

    if (!save_dir.empty()) {
        std::string contact_scaling_string = std::to_string(contact_scaling);
        std::string tmax_string            = std::to_string(tmax);
        std::string dt_string              = std::to_string(simulation_parameter["dt"]);

        std::string filename_ide =
            save_dir + "changepoint_ide_" + contact_scaling_string.substr(0, contact_scaling_string.find(".") + 2) +
            "_" + tmax_string.substr(0, tmax_string.find(".")) + "_" + dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ide_flows = filename_ide + "_flows.h5";
        mio::IOResult<void> save_result_status_f =
            mio::save_result({sim.get_transitions()}, {0}, 1, filename_ide_flows);

        std::string filename_ide_compartments = filename_ide + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({sim.get_result()}, {0}, 1, filename_ide_compartments);
    }

    // Return vector with initial compartments.
    return mio::success(sim.get_result());
}

/**
* @brief Function that simulates from time 0 until tmax using an ODE model where we apply a contact scaling after
* two days.   
*
* @param[in] init_compartments Vector containing initial values for the compartments. 
* @param[in] contact_scaling Factor that is applied to contact matrix after two days. 
* @param[in] tmax Time up to which we simulate. 
* @param[in] save_dir Directory where simulation results will be stored. 
* @returns Any io errors that happen.
*/
mio::IOResult<void> simulate_ode_model(Vector init_compartments, ScalarType contact_scaling, ScalarType tmax,
                                       std::string save_dir = "")
{
    // Use ODE FlowModel.
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
        1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
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

    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() = scale_contact_matrix(contact_scaling);

    model_ode.check_constraints();

    // Set integrator and fix step size.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    integrator->set_dt_min(simulation_parameter["dt"]);
    integrator->set_dt_max(simulation_parameter["dt"]);

    // Simulate.
    std::vector<mio::TimeSeries<ScalarType>> results_ode = mio::osecir::simulate_flows<ScalarType>(
        simulation_parameter["t0"], tmax, simulation_parameter["dt"], model_ode, integrator);

    // Save results.
    if (!save_dir.empty()) {
        std::string contact_scaling_string = std::to_string(contact_scaling);
        std::string tmax_string            = std::to_string(tmax);
        std::string dt_string              = std::to_string(simulation_parameter["dt"]);

        std::string filename_ode =
            save_dir + "changepoint_ode_" + contact_scaling_string.substr(0, contact_scaling_string.find(".") + 2) +
            "_" + tmax_string.substr(0, tmax_string.find(".")) + "_" + dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ode_flows           = filename_ode + "_flows.h5";
        mio::IOResult<void> save_result_status_f = mio::save_result({results_ode[1]}, {0}, 1, filename_ode_flows);

        std::string filename_ode_compartments = filename_ode + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({results_ode[0]}, {0}, 1, filename_ode_compartments);
    }

    return mio::success();
}

int main()
{
    // Paths are valid if file is executed e.g. in memilio/build/bin.
    std::string save_dir = "../../data/simulation_results/changepoints/";
    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    // Changepoint scenario with halving of contacts after two days.
    ScalarType contact_scaling = 0.5;
    ScalarType tmax            = 12;

    auto result_ide = simulate_ide_model(contact_scaling, tmax, save_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    // Use compartments at time 0 from IDE simulation as initial values for ODE model to make results comparable.
    Vector compartments = result_ide.value().get_value(0);

    auto result_ode = simulate_ode_model(compartments, contact_scaling, tmax, save_dir);
    if (!result_ode) {
        printf("%s\n", result_ode.error().formatted_message().c_str());
        return -1;
    }

    // Changepoint scenario with doubling of contacts after two days.
    contact_scaling = 2.;
    tmax            = 12;

    result_ide = simulate_ide_model(contact_scaling, tmax, save_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    // Use compartments at time 0 from IDE simulation as initial values for ODE model to make results comparable.
    compartments = result_ide.value().get_value(0);

    result_ode = simulate_ode_model(compartments, contact_scaling, tmax, save_dir);
    if (!result_ode) {
        printf("%s\n", result_ode.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
