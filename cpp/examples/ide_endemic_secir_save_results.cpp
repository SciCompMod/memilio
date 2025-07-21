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

std::map<std::string, ScalarType> simulation_parameter = {{"dt", 1.},
                                                          {"t0", 0.},
                                                          {"Susceptibles", 100000.},
                                                          {"Exposed", 0},
                                                          {"InfectedNoSymptoms", 10.},
                                                          {"InfectedSymptoms", 20.},
                                                          {"InfectedSevere", 0.},
                                                          {"InfectedCritical", 0.},
                                                          {"Recovered", 0.},
                                                          {"Dead", 0.},
                                                          {"TransmissionProbabilityOnContact", 0.1},
                                                          {"RelativeTransmissionNoSymptoms", 0.5},
                                                          {"RiskOfInfectionFromSymptomatic", 0.5},
                                                          {"InfectedSymptomsPerInfectedNoSymptoms", 0.5},
                                                          {"SeverePerInfectedSymptoms", 0.1},
                                                          {"CriticalPerSevere", 0.3},
                                                          {"DeathsPerCritical", 0.2},
                                                          {"BirthRate1", 4e-3},
                                                          {"DeathRate1", 3e-3},
                                                          {"BirthRate2", 3e-4},
                                                          {"DeathRate2", 3e-4},
                                                          {"BirthRate3", 3e-4},
                                                          {"DeathRate3", 4e-4},
                                                          {"Contacts", 10.}};

mio::IOResult<void> simulate_endidemodel(ScalarType tmax, std::string save_dir = "")
{
    int num_states      = static_cast<int>(mio::endisecir::InfectionState::Count);
    int num_transitions = static_cast<int>(mio::endisecir::InfectionTransition::Count);

    //Set the inital values for the compartments.
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

    // We are going to simulate thre times with three different Birth and Death Rates,
    // but all other parameters stay the same.

    // Model 1, with BirthRate > DeathRate:
    mio::TimeSeries<ScalarType> init_copy1(init);
    mio::endisecir::Model model_endide1(std::move(init_copy1));

    // Set distributions.
    // Set TransitionDistributions.

    mio::SmootherCosine smoothcos(8.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    // Uncomment for Lognorm.
    // mio::ConstantFunction initialfunc(0);
    // mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution((int)mio::endisecir::InfectionTransition::Count,
    //                                                                 delaydistributioninit);
    // // ExposedToInfectedNoSymptoms
    // mio::LognormSurvivalFunction survivalExposedToInfectedNoSymptoms(0.3, 0, 4.2);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
    //     survivalExposedToInfectedNoSymptoms);
    // // InfectedNoSymptomsToInfectedSymptoms
    // mio::LognormSurvivalFunction survivalInfectedNoSymptomsToInfectedSymptoms(0.7, 0, 0.8);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
    //     .set_state_age_function(survivalInfectedNoSymptomsToInfectedSymptoms);
    // // InfectedNoSymptomsToRecovered
    // mio::LognormSurvivalFunction survivalInfectedNoSymptomsToRecovered(0.2, 0, 7.7);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedNoSymptomsToRecovered]
    //     .set_state_age_function(survivalInfectedNoSymptomsToRecovered);
    // // InfectedSymptomsToInfectedSevere
    // mio::LognormSurvivalFunction survivalInfectedSymptomsToInfectedSevere(0.7, 0, 5.3);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
    //     .set_state_age_function(survivalInfectedSymptomsToInfectedSevere);
    // // InfectedSymptomsToRecovered
    // mio::LognormSurvivalFunction survivalInfectedSymptomsToRecovered(0.2, 0, 7.8);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedSymptomsToRecovered].set_state_age_function(
    //     survivalInfectedSymptomsToRecovered);
    // // InfectedSevereToInfectedCritical
    // mio::LognormSurvivalFunction survivalInfectedSevereToInfectedCritical(1.0, 0, 0.9);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedSevereToInfectedCritical]
    //     .set_state_age_function(survivalInfectedSevereToInfectedCritical);
    // // InfectedSevereToRecovered
    // mio::LognormSurvivalFunction survivalInfectedSevereToRecovered(0.3, 0, 17.1);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedSevereToRecovered].set_state_age_function(
    //     survivalInfectedSevereToRecovered);
    // // InfectedCriticalToDead
    // mio::LognormSurvivalFunction survivalInfectedCriticalToDead(0.4, 0, 9.8);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedCriticalToDead].set_state_age_function(
    //     survivalInfectedCriticalToDead);
    // // InfectedCriticalToRecovered
    // mio::LognormSurvivalFunction survivalInfectedCriticalToRecovered(0.3, 0, 17.1);
    // vec_delaydistribution[(int)mio::endisecir::InfectionTransition::InfectedCriticalToRecovered].set_state_age_function(
    //     survivalInfectedCriticalToRecovered);

    model_endide1.parameters.set<mio::endisecir::TransitionDistributions>(vec_delaydistribution);

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

    mio::ConstantFunction constant(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper constant_prob(constant);

    model_endide1.parameters.get<mio::endisecir::TransmissionProbabilityOnContact>() = constant_prob;

    mio::ExponentialSurvivalFunction exponential(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    mio::StateAgeFunctionWrapper exponential_prob(exponential);

    model_endide1.parameters.set<mio::endisecir::TransmissionProbabilityOnContact>(exponential_prob);
    exponential_prob.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_endide1.parameters.set<mio::endisecir::RiskOfInfectionFromSymptomatic>(exponential_prob);

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

    model_endide2.parameters.get<mio::endisecir::TransmissionProbabilityOnContact>() = constant_prob;

    exponential_prob.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_endide2.parameters.set<mio::endisecir::TransmissionProbabilityOnContact>(exponential_prob);
    exponential_prob.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_endide2.parameters.set<mio::endisecir::RiskOfInfectionFromSymptomatic>(exponential_prob);

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

    model_endide3.parameters.get<mio::endisecir::TransmissionProbabilityOnContact>() = constant_prob;

    exponential_prob.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_endide3.parameters.set<mio::endisecir::TransmissionProbabilityOnContact>(exponential_prob);
    exponential_prob.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_endide3.parameters.set<mio::endisecir::RiskOfInfectionFromSymptomatic>(exponential_prob);

    model_endide3.parameters.set<mio::endisecir::NaturalBirthRate>(simulation_parameter["BirthRate3"]);
    model_endide3.parameters.set<mio::endisecir::NaturalDeathRate>(simulation_parameter["DeathRate3"]);

    model_endide3.set_tol_for_support_max(1e-6);
    model_endide3.check_constraints();

    // Simulate.
    mio::endisecir::Simulation sim1(model_endide1, simulation_parameter["dt"]);
    sim1.advance(tmax);
    // mio::endisecir::Simulation sim2(model_endide2, simulation_parameter["dt"]);
    // sim2.advance(tmax);
    // mio::endisecir::Simulation sim3(model_endide3, simulation_parameter["dt"]);
    // sim3.advance(tmax);
    sim1.get_forceofinfections().print_table({"FoI"}, 16, 8);

    if (!save_dir.empty()) {
        std::string tmax_string = std::to_string(tmax);
        std::string dt_string   = std::to_string(simulation_parameter["dt"]);

        std::string filename_ide = save_dir + "analysis_endide_" + tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        //Save total population for different birth and death rates.
        std::string filename_ide_totalpopulation1 = filename_ide + "_totalpopulation1.h5";
        mio::IOResult<void> save_result_status_tp1 =
            mio::save_result({sim1.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation1);

        // std::string filename_ide_totalpopulation2 = filename_ide + "_totalpopulation2.h5";
        // mio::IOResult<void> save_result_status_tp2 =
        //     mio::save_result({sim2.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation2);

        // std::string filename_ide_totalpopulation3 = filename_ide + "_totalpopulation3.h5";
        // mio::IOResult<void> save_result_status_tp3 =
        //     mio::save_result({sim3.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation3);

        //Save compartments for different birth and death rates.
        std::string filename_ide_compartments1 = filename_ide + "_compartments1.h5";
        mio::IOResult<void> save_result_status_c1 =
            mio::save_result({sim1.get_compartments()}, {0}, 1, filename_ide_compartments1);

        // std::string filename_ide_compartments2 = filename_ide + "_compartments2.h5";
        // mio::IOResult<void> save_result_status_c2 =
        //     mio::save_result({sim2.get_compartments()}, {0}, 1, filename_ide_compartments2);

        // std::string filename_ide_compartments3 = filename_ide + "_compartments3.h5";
        // mio::IOResult<void> save_result_status_c3 =
        //     mio::save_result({sim3.get_compartments()}, {0}, 1, filename_ide_compartments3);

        //Save normalized compartments for only one birth and death rate.
        std::string filename_ide_normcompartments = filename_ide + "_normcompartments.h5";
        mio::IOResult<void> save_result_status_nc =
            mio::save_result({sim1.get_normalizedcompartments()}, {0}, 1, filename_ide_normcompartments);

        //Save the equilibrium of the normalized compartments.
        // std::string filename_ide_normcompartmentsequi = filename_ide + "_equinormcompartments.h5";
        // mio::IOResult<void> save_result_status_nce =
        //     mio::save_result({sim1.get_equilibriumcompartments()}, {0}, 1, filename_ide_normcompartmentsequi);

        //Save the force of infection for only one birth and death rate.
        std::string filename_ide_forceofinfection = filename_ide + "_forceofinfection.h5";
        mio::IOResult<void> save_result_status_foi =
            mio::save_result({sim1.get_forceofinfections()}, {0}, 1, filename_ide_forceofinfection);

        //Save the equilibrium of the force of infection.
        // std::string filename_ide_forceofinfectionequi = filename_ide + "_equiforceofinfection.h5";
        // mio::IOResult<void> save_result_status_foie =
        //     mio::save_result({sim1.get_equilibrium_forceofinfection()}, {0}, 1, filename_ide_forceofinfectionequi);

        // if (!save_result_status_tp1 || !save_result_status_tp2 || !save_result_status_tp3 || !save_result_status_c1 ||
        //     !save_result_status_c2 || !save_result_status_c3 || !save_result_status_nc || !save_result_status_foi ||
        //     !save_result_status_nce || !save_result_status_foie) {
        //     return mio::failure(mio::StatusCode::UnknownError, "Error while saving results.");
        // }
        if (!save_result_status_tp1 || !save_result_status_c1 || !save_result_status_nc || !save_result_status_foi) {
            return mio::failure(mio::StatusCode::UnknownError, "Error while saving results.");
        }
    }
    // Print the reproduction numbers for different birth and death rates.
    std::cout << "The reproduction number Rc for Birth rate > Death rate is  " << sim1.get_reproductionnumber_c()
              << "\n";
    // std::cout << "The reproduction number Rc for Birth rate = Death rate is  " << sim2.get_reproductionnumber_c()
    //           << "\n";
    // std::cout << "The reproduction number Rc for Birth rate < Death rate is  " << sim3.get_reproductionnumber_c()
    //           << "\n";

    // Return results of the simulation.
    return mio::success();
}

int main()
{
    std::string result_dir = "/localdata1/trit_ha/code/memilio-1/PythonPlotsEndIDE/simulation_results/";

    // Define tmax.
    ScalarType tmax = 100;

    auto result_ide = simulate_endidemodel(tmax, result_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}