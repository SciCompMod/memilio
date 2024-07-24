/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler
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
#include "ide_secir/parameters_io.h"
#include "ide_secir/simulation.h"

#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"

#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "ode_secir/infection_state.h"
#include <string>
#include <map>
#include <iostream>

using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"t0", 0.},
                                                          {"dt_flows", 0.1},
                                                          {"total_population", 83155031.},
                                                          {"total_confirmed_cases", 0.}, //341223.
                                                          {"deaths", 9710.},
                                                          {"TimeExposed", 4.5},
                                                          {"TimeInfectedNoSymptoms", 2.52761690},
                                                          {"TimeInfectedSymptoms", 7.88989994},
                                                          {"TimeInfectedSevere", 15.22527840},
                                                          {"TimeInfectedCritical", 16.49289020},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"RecoveredPerInfectedNoSymptoms", 0.206901},
                                                          {"SeverePerInfectedSymptoms", 0.0786429},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.217177},
                                                          {"cont_freq", 5.5}};

std::vector<mio::TimeSeries<ScalarType>> simulate_ide_model(std::string date, ScalarType tmax,
                                                            std::string save_dir = "")
{
    // Initialize model.
    mio::isecir::Model model_ide(mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionTransition::Count),
                                 simulation_parameter["total_population"], simulation_parameter["deaths"],
                                 simulation_parameter["total_confirmed_cases"]);

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

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, simulation_parameter["cont_freq"]));
    model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

    model_ide.set_tol_for_support_max(1e-6);

    const fs::path& data_dir = "../../data";

    std::string path_rki = mio::path_join((data_dir / "pydata" / "Germany").string(), "cases_all_germany.json");

    mio::Date parsed_date = mio::parse_date(date).value();
    mio::IOResult<void> init_flows =
        set_initial_flows(model_ide, simulation_parameter["dt_flows"], path_rki, parsed_date, 1.9);

    model_ide.check_constraints(simulation_parameter["dt_flows"]);

    std::cout << "Global support max: " << model_ide.get_global_support_max(simulation_parameter["dt_flows"])
              << std::endl;

    // Simulate.
    mio::isecir::Simulation sim(model_ide, simulation_parameter["dt_flows"]);
    sim.get_transitions().print_table();
    sim.advance(tmax);

    if (!save_dir.empty()) {
        std::string tmax_string  = std::to_string(tmax);
        std::string dt_string    = std::to_string(simulation_parameter["dt_flows"]);
        std::string filename_ide = save_dir + "ide_" + date + "_" + tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ide_flows = filename_ide + "_flows.h5";
        mio::IOResult<void> save_result_status_f =
            mio::save_result({sim.get_transitions()}, {0}, 1, filename_ide_flows);
        std::string filename_ide_compartments = filename_ide + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({sim.get_result()}, {0}, 1, filename_ide_compartments);
    }

    // Return vector with simulation results.
    return {sim.get_result(), sim.get_transitions()};
}

mio::IOResult<void> simulate_ode_model(std::string date, Vector init_compartments, ScalarType tmax,
                                       std::string save_dir = "")
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

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, simulation_parameter["cont_freq"]));
    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    model_ode.check_constraints();

    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    integrator->set_dt_min(simulation_parameter["dt_flows"]);
    integrator->set_dt_max(simulation_parameter["dt_flows"]);

    std::vector<mio::TimeSeries<ScalarType>> results_ode = mio::osecir::simulate_flows<ScalarType>(
        simulation_parameter["t0"], tmax, simulation_parameter["dt_flows"], model_ode, integrator);

    if (!save_dir.empty()) {
        std::string tmax_string  = std::to_string(tmax);
        std::string dt_string    = std::to_string(simulation_parameter["dt_flows"]);
        std::string filename_ode = save_dir + "ode_" + date + "_" + tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ode_flows           = filename_ode + "_flows.h5";
        mio::IOResult<void> save_result_status_f = mio::save_result({results_ode[1]}, {0}, 1, filename_ode_flows);
        std::string filename_ode_compartments    = filename_ode + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({results_ode[0]}, {0}, 1, filename_ode_compartments);
    }

    return mio::success();
}

int main()
{
    // Paths are valid if file is executed eg in memilio/build/bin.
    std::string save_dir = "../../results/real/";
    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    std::string date = "2020-10-01";
    ScalarType tmax  = 30;

    std::vector<mio::TimeSeries<ScalarType>> result_ide = simulate_ide_model(date, tmax, save_dir);

    // result_ide[1].print_table();

    Vector init_compartments = result_ide[0].get_value(0);

    auto result_ode = simulate_ode_model(date, init_compartments, tmax, save_dir);
    if (!result_ode) {
        printf("%s\n", result_ode.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}