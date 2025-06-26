#include "memilio/io/result_io.h"
#include "memilio/utils/time_series.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/floating_point.h"

#include "ide_endemic_secir/infection_state.h"
#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/simulation.h"
#include "ide_endemic_secir/parameters.h"

#include "boost/filesystem.hpp"
#include <iostream>
#include <string>
#include <utility>

std::map<std::string, ScalarType> simulation_parameter = {{"dt", 0.1},
                                                          {"t0", 0.},
                                                          {"Susceptibles", 10000.},
                                                          {"Exposed", 0},
                                                          {"InfectedNoSymptoms", 0.},
                                                          {"InfectedSymptoms", 0.},
                                                          {"InfectedSevere", 0.},
                                                          {"InfectedCritical", 0.},
                                                          {"Recovered", 0.},
                                                          {"Dead", 0.},
                                                          {"TransmissionProbabilityOnContact", 0.5},
                                                          {"RelativeTransmissionNoSymptoms", 0.5},
                                                          {"RiskOfInfectionFromSymptomatic", 0.5},
                                                          {"InfectedSymptomsPerInfectedNoSymptoms", 0.5},
                                                          {"SeverePerInfectedSymptoms", 0.5},
                                                          {"CriticalPerSevere", 0.5},
                                                          {"DeathsPerCritical", 0.1},
                                                          {"BirthRate1", 4e-5},
                                                          {"DeathRate1", 3e-5},
                                                          {"BirthRate2", 3e-5},
                                                          {"DeathRate2", 3e-5},
                                                          {"BirthRate3", 3e-5},
                                                          {"DeathRate3", 4e-5},
                                                          {"Contacts", 10.}};

mio::IOResult<void> simulate_endidemodel(ScalarType tmax, std::string save_dir = "")
{
    int num_states      = static_cast<int>(mio::endisecir::InfectionState::Count);
    int num_transitions = static_cast<int>(mio::endisecir::InfectionTransition::Count);

    mio::TimeSeries<ScalarType> init(num_states);

    mio::TimeSeries<ScalarType>::Vector vec_init(num_states);

    vec_init[static_cast<int>(mio::endisecir::InfectionState::Susceptible)] = simulation_parameter["Susceptibles"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Exposed)]     = simulation_parameter["Exposed"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedNoSymptoms)] =
        simulation_parameter["InfectedNoSymptoms"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedSymptoms)] =
        simulation_parameter["InfectedSymptoms"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedSevere)] = simulation_parameter["InfectedSevere"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedCritical)] =
        simulation_parameter["InfectedCritical"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Recovered)] = simulation_parameter["Recovered"];
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Dead)]      = simulation_parameter["Dead"];

    init.add_time_point(0, vec_init);

    // Model 1, with BirthRate > DeathRate:
    mio::TimeSeries<ScalarType> init_copy1(init);
    mio::endisecir::Model model_endide1(std::move(init_copy1));

    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    model_endide1.parameters.get<mio::endisecir::TransitionDistributions>() = vec_delaydistribution;

    // Set other parameters.
    std::vector<ScalarType> vec_prob((int)mio::endisecir::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
        simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
        1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSymptomsToRecovered)] =
        1 - simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
        simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSevereToRecovered)] =
        1 - simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedCriticalToDead)] =
        simulation_parameter["DeathsPerCritical"];
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedCriticalToRecovered)] =
        1 - simulation_parameter["DeathsPerCritical"];

    model_endide1.parameters.set<mio::endisecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, simulation_parameter["Contacts"]));
    model_endide1.parameters.get<mio::endisecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
    model_endide1.parameters.set<mio::endisecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_endide1.parameters.set<mio::endisecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_endide1.parameters.set<mio::endisecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

    model_endide1.parameters.set<mio::endisecir::NaturalBirthRate>(simulation_parameter["BirthRate1"]);
    model_endide1.parameters.set<mio::endisecir::NaturalDeathRate>(simulation_parameter["DeathRate1"]);

    model_endide1.set_tol_for_support_max(1e-6);
    model_endide1.check_constraints();

    // Model 2, with BirthRate = DeathRate:
    mio::TimeSeries<ScalarType> init_copy2(init);
    mio::endisecir::Model model_endide2(std::move(init_copy2));

    model_endide2.parameters.get<mio::endisecir::TransitionDistributions>() = vec_delaydistribution;

    // Set other parameters.
    model_endide2.parameters.set<mio::endisecir::TransitionProbabilities>(vec_prob);

    model_endide2.parameters.get<mio::endisecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc2(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide2(constfunc);
    model_endide2.parameters.set<mio::endisecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide2);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_endide2.parameters.set<mio::endisecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide2);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_endide2.parameters.set<mio::endisecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide2);

    model_endide2.parameters.set<mio::endisecir::NaturalBirthRate>(simulation_parameter["BirthRate2"]);
    model_endide2.parameters.set<mio::endisecir::NaturalDeathRate>(simulation_parameter["DeathRate2"]);

    model_endide2.set_tol_for_support_max(1e-6);
    model_endide2.check_constraints();

    // Model 3, with BirthRate < DeathRate:
    mio::TimeSeries<ScalarType> init_copy3(init);
    mio::endisecir::Model model_endide3(std::move(init_copy3));

    model_endide3.parameters.get<mio::endisecir::TransitionDistributions>() = vec_delaydistribution;

    // Set other parameters.
    model_endide3.parameters.set<mio::endisecir::TransitionProbabilities>(vec_prob);

    model_endide3.parameters.get<mio::endisecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constfunc3(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide3(constfunc);
    model_endide3.parameters.set<mio::endisecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide3);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_endide3.parameters.set<mio::endisecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide3);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_endide3.parameters.set<mio::endisecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide3);

    model_endide3.parameters.set<mio::endisecir::NaturalBirthRate>(simulation_parameter["BirthRate3"]);
    model_endide3.parameters.set<mio::endisecir::NaturalDeathRate>(simulation_parameter["DeathRate3"]);

    model_endide3.set_tol_for_support_max(1e-6);
    model_endide3.check_constraints();

    // Simulate.
    mio::endisecir::Simulation sim1(model_endide1, simulation_parameter["dt"]);
    sim1.advance(tmax);
    mio::endisecir::Simulation sim2(model_endide2, simulation_parameter["dt"]);
    sim2.advance(tmax);
    mio::endisecir::Simulation sim3(model_endide3, simulation_parameter["dt"]);
    sim3.advance(tmax);

    if (!save_dir.empty()) {
        std::string tmax_string = std::to_string(tmax);
        std::string dt_string   = std::to_string(simulation_parameter["dt"]);

        std::string filename_ide = save_dir + "analysis_endide_" + tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        //Save total population for different birth and death rates
        std::string filename_ide_totalpopulation1 = filename_ide + "_totalpopulation1.h5";
        mio::IOResult<void> save_result_status_tp1 =
            mio::save_result({sim1.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation1);

        std::string filename_ide_totalpopulation2 = filename_ide + "_totalpopulation2.h5";
        mio::IOResult<void> save_result_status_tp2 =
            mio::save_result({sim2.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation2);

        std::string filename_ide_totalpopulation3 = filename_ide + "_totalpopulation3.h5";
        mio::IOResult<void> save_result_status_tp3 =
            mio::save_result({sim3.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation3);

        //Save compartments for different birth and death rates
        std::string filename_ide_compartments1 = filename_ide + "_compartments1.h5";
        mio::IOResult<void> save_result_status_c1 =
            mio::save_result({sim1.get_compartments()}, {0}, 1, filename_ide_compartments1);

        std::string filename_ide_compartments2 = filename_ide + "_compartments2.h5";
        mio::IOResult<void> save_result_status_c2 =
            mio::save_result({sim2.get_compartments()}, {0}, 1, filename_ide_compartments2);

        std::string filename_ide_compartments3 = filename_ide + "_compartments3.h5";
        mio::IOResult<void> save_result_status_c3 =
            mio::save_result({sim3.get_compartments()}, {0}, 1, filename_ide_compartments3);

        if (!save_result_status_tp1 || !save_result_status_tp2 || !save_result_status_tp3 || !save_result_status_c1 ||
            !save_result_status_c2 || !save_result_status_c3) {
            return mio::failure(mio::StatusCode::UnknownError, "Error while saving results.");
        }
    }

    // Return results (i.e. total populations) of the simulation.
    return mio::success();
}

int main()
{
    std::string result_dir = "/localdata1/trit_ha/code/memilio-1/PythonPlotsEndIDE/simulation_results/";

    // Define tmax for both scenarios.
    ScalarType tmax = 500;

    auto result_ide = simulate_endidemodel(tmax, result_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}