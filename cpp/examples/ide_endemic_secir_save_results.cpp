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

mio::IOResult<void> simulate_endidemodel(ScalarType tmax, std::string save_dir = "")
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType dt = 1.0;

    int num_states      = static_cast<int>(mio::endisecir::InfectionState::Count);
    int num_transitions = static_cast<int>(mio::endisecir::InfectionTransition::Count);

    // Create TimeSeries with num_states elements where states needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_states);

    Vec vec_init(num_states);

    vec_init[static_cast<int>(mio::endisecir::InfectionState::Susceptible)]        = 100000.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Exposed)]            = 0.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedNoSymptoms)] = 10.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedSymptoms)]   = 20.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedSevere)]     = 0.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedCritical)]   = 0.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Recovered)]          = 0.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Dead)]               = 0.;

    init.add_time_point(0, vec_init);

    mio::endisecir::CompParameters computed_parameters(std::move(init));

    //Set working parameters

    // mio::ExponentialSurvivalFunction exp(3.0);
    // mio::StateAgeFunctionWrapper delaydistribution(exp);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    mio::SmootherCosine smoothcos(6.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    //Uncomment for Lognorm.
    // mio::ConstantFunction initialfunc(0);
    // mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistributioninit);
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

    computed_parameters.parameters.get<mio::endisecir::TransitionDistributions>() = vec_delaydistribution;

    std::vector<ScalarType> vec_prob((int)mio::endisecir::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] = 0.8;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedNoSymptomsToRecovered)]        = 1 - 0.8;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]     = 0.1;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSymptomsToRecovered)]          = 1 - 0.1;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSevereToInfectedCritical)]     = 0.2;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSevereToRecovered)]            = 1 - 0.2;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedCriticalToDead)]               = 0.4;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedCriticalToRecovered)]          = 1 - 0.4;

    //Different values:
    // std::vector<ScalarType> vec_prob((int)mio::endisecir::InfectionTransition::Count, 1.);
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] = 0.5;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedNoSymptomsToRecovered)]        = 1 - 0.5;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]     = 0.5;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSymptomsToRecovered)]          = 1 - 0.5;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSevereToInfectedCritical)]     = 0.5;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedSevereToRecovered)]            = 1 - 0.5;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedCriticalToDead)]               = 0.8;
    // vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::InfectedCriticalToRecovered)]          = 1 - 0.8;

    computed_parameters.parameters.set<mio::endisecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    computed_parameters.parameters.get<mio::endisecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constant(0.1);
    mio::StateAgeFunctionWrapper constant_prob(constant);

    computed_parameters.parameters.get<mio::endisecir::TransmissionProbabilityOnContact>() = constant_prob;

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper exponential_prob(exponential);

    computed_parameters.parameters.get<mio::endisecir::RelativeTransmissionNoSymptoms>() = exponential_prob;
    computed_parameters.parameters.get<mio::endisecir::RiskOfInfectionFromSymptomatic>() = exponential_prob;

    computed_parameters.parameters.set<mio::endisecir::NaturalBirthRate>(4e-3);
    computed_parameters.parameters.set<mio::endisecir::NaturalDeathRate>(3e-3);

    //computed_parameters.set_tol_for_support_max(1e-6);

    mio::endisecir::Model model(computed_parameters);

    mio::endisecir::NormModel normmodel(computed_parameters);

    // start the simulation.
    mio::endisecir::Simulation sim(computed_parameters, model, normmodel, dt);
    sim.advance(tmax);

    if (!save_dir.empty()) {
        std::string tmax_string = std::to_string(tmax);
        std::string dt_string   = std::to_string(dt);

        std::string filename_ide = save_dir + "analysis_endide_" + tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        //Save files of Model.
        //Save total population.
        std::string filename_ide_totalpopulation1 = filename_ide + "_totalpopulation1.h5";
        mio::IOResult<void> save_result_status_tp1 =
            mio::save_result({sim.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation1);
        //Save derivative of total population.
        std::string filename_ide_dertotalpopulation1 = filename_ide + "_dertotalpopulation1.h5";
        mio::IOResult<void> save_result_status_dtp1 =
            mio::save_result({sim.get_totalpopulations_derivative()}, {0}, 1, filename_ide_dertotalpopulation1);

        //Save compartments.
        std::string filename_ide_compartments1 = filename_ide + "_compartments1.h5";
        mio::IOResult<void> save_result_status_c1 =
            mio::save_result({sim.get_compartments()}, {0}, 1, filename_ide_compartments1);

        //Save normalized compartments.
        std::string filename_ide_normcompartments = filename_ide + "_normcompartments.h5";
        mio::IOResult<void> save_result_status_nc =
            mio::save_result({sim.get_normalizedcompartments()}, {0}, 1, filename_ide_normcompartments);

        //Save the force of infection.
        std::string filename_ide_forceofinfection = filename_ide + "_forceofinfection.h5";
        mio::IOResult<void> save_result_status_foi =
            mio::save_result({sim.get_forceofinfections()}, {0}, 1, filename_ide_forceofinfection);

        //Save files of NormModel.
        //Save compartments.
        std::string filename_ide_normmod_compartments = filename_ide + "_normmod_compartments.h5";
        mio::IOResult<void> save_result_status_nmc1 =
            mio::save_result({sim.get_normmodel_compartments()}, {0}, 1, filename_ide_normmod_compartments);

        //Save the force of infection.
        std::string filename_ide_normmod_forceofinfection = filename_ide + "_normmod_forceofinfection.h5";
        mio::IOResult<void> save_result_status_nmfoi =
            mio::save_result({sim.get_normmodel_forceofinfections()}, {0}, 1, filename_ide_normmod_forceofinfection);

        //Safe difference between normalized compartments.
        std::string filename_ide_difference_normcomp = filename_ide + "_difference_normcomp.h5";
        mio::IOResult<void> save_result_status_dnm =
            mio::save_result({sim.get_difference_normalizationcomp()}, {0}, 1, filename_ide_difference_normcomp);

        //Safe difference between Force of Infections.
        std::string filename_ide_difference_normFoI = filename_ide + "_difference_normfoi.h5";
        mio::IOResult<void> save_result_status_dnf =
            mio::save_result({sim.get_difference_normalizationFoi()}, {0}, 1, filename_ide_difference_normFoI);

        if (!save_result_status_tp1 || !save_result_status_c1 || !save_result_status_nc || !save_result_status_foi ||
            !save_result_status_nmc1 || !save_result_status_nmfoi || !save_result_status_dnm ||
            !save_result_status_dnf || !save_result_status_dtp1) {
            return mio::failure(mio::StatusCode::UnknownError, "Error while saving results.");
        }
    }
    // Print the reproduction number.
    std::cout << "The reproduction number Rc for Birth rate > Death rate is  " << sim.get_reproductionnumber_c()
              << "\n";

    return mio::success();
}

int main()
{
    std::string result_dir = "/localdata1/trit_ha/code/memilio-1/PythonPlotsEndIDESECIR/simulation_results/";

    // Define tmax.
    ScalarType tmax = 2000;

    auto result_ide = simulate_endidemodel(tmax, result_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}