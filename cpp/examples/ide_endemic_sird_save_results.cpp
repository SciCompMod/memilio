#include "memilio/io/result_io.h"
#include "memilio/utils/time_series.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/floating_point.h"

#include "ide_endemic_sird/infection_state.h"
#include "ide_endemic_sird/model.h"
#include "ide_endemic_sird/simulation.h"
#include "ide_endemic_sird/parameters.h"

#include "boost/filesystem.hpp"
#include <iostream>
#include <string>
#include <utility>

mio::IOResult<void> simulate_endidemodel(ScalarType tmax, std::string save_dir = "")
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType dt = 1.0;

    int num_states      = static_cast<int>(mio::endisird::InfectionState::Count);
    int num_transitions = static_cast<int>(mio::endisird::InfectionTransition::Count);

    // Create TimeSeries with num_states elements where states needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_states);

    Vec vec_init(num_states);

    vec_init[static_cast<int>(mio::endisird::InfectionState::Susceptible)] = 100000.;
    vec_init[static_cast<int>(mio::endisird::InfectionState::Infected)]    = 30.;
    vec_init[static_cast<int>(mio::endisird::InfectionState::Recovered)]   = 0.;
    vec_init[static_cast<int>(mio::endisird::InfectionState::Dead)]        = 0.;

    init.add_time_point(0, vec_init);

    mio::endisird::CompParameters computed_parameters(std::move(init));

    //Set working parameters

    // mio::ExponentialSurvivalFunction exp(3.0);
    // mio::StateAgeFunctionWrapper delaydistribution(exp);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    mio::SmootherCosine smoothcos(12.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    // mio::ConstantFunction initialfunc(0);
    // mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistributioninit);
    // // InfectedToDead
    // mio::SmootherCosine survivalExposedToInfectedNoSymptoms(14.0);
    // vec_delaydistribution[(int)mio::endisird::InfectionTransition::InfectedToDead].set_state_age_function(
    //     survivalExposedToInfectedNoSymptoms);
    // // InfectedToRecovered
    // mio::SmootherCosine survivalInfectedNoSymptomsToInfectedSymptoms(8.0);
    // vec_delaydistribution[(int)mio::endisird::InfectionTransition::InfectedToRecovered].set_state_age_function(
    //     survivalInfectedNoSymptomsToInfectedSymptoms);

    computed_parameters.parameters.get<mio::endisird::TransitionDistributions>() = vec_delaydistribution;

    std::vector<ScalarType> vec_prob((int)mio::endisird::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::endisird::InfectionTransition::InfectedToDead)]      = 0.1;
    vec_prob[Eigen::Index(mio::endisird::InfectionTransition::InfectedToRecovered)] = 1 - 0.1;

    computed_parameters.parameters.set<mio::endisird::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    computed_parameters.parameters.get<mio::endisird::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constant(0.1);
    mio::StateAgeFunctionWrapper constant_prob(constant);

    computed_parameters.parameters.get<mio::endisird::TransmissionProbabilityOnContact>() = constant_prob;

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper exponential_prob(exponential);

    computed_parameters.parameters.get<mio::endisird::RiskOfInfectionFromSymptomatic>() = exponential_prob;

    computed_parameters.parameters.set<mio::endisird::NaturalBirthRate>(1e-1);
    computed_parameters.parameters.set<mio::endisird::NaturalDeathRate>(2e-2);

    //computed_parameters.set_tol_for_support_max(1e-6);

    mio::endisird::Model model(computed_parameters);

    mio::endisird::NormModel normmodel(computed_parameters);

    // start the simulation.
    mio::endisird::Simulation sim(computed_parameters, model, normmodel, dt);
    sim.advance(tmax);

    if (!save_dir.empty()) {
        std::string tmax_string = std::to_string(tmax);
        std::string dt_string   = std::to_string(dt);

        std::string filename_ide = save_dir + "analysis_endide_SIRD_" + tmax_string.substr(0, tmax_string.find(".")) +
                                   "_" + dt_string.substr(0, dt_string.find(".") + 5);

        //Save files of Model.
        //Save total population.
        std::string filename_ide_totalpopulation1 = filename_ide + "_totalpopulation.h5";
        mio::IOResult<void> save_result_status_tp1 =
            mio::save_result({sim.get_totalpopulations()}, {0}, 1, filename_ide_totalpopulation1);
        //Save derivative of total population.
        std::string filename_ide_dertotalpopulation1 = filename_ide + "_dertotalpopulation.h5";
        mio::IOResult<void> save_result_status_dtp1 =
            mio::save_result({sim.get_totalpopulations_derivative()}, {0}, 1, filename_ide_dertotalpopulation1);

        //Save compartments.
        std::string filename_ide_compartments1 = filename_ide + "_compartments.h5";
        mio::IOResult<void> save_result_status_c1 =
            mio::save_result({sim.get_compartments()}, {0}, 1, filename_ide_compartments1);

        //Save normalized compartments.
        std::string filename_ide_normcompartments = filename_ide + "_normcompartments.h5";
        mio::IOResult<void> save_result_status_nc =
            mio::save_result({sim.get_normalizedcompartments()}, {0}, 1, filename_ide_normcompartments);

        //Save the force of infection.
        std::string filename_ide_forceofinfection = filename_ide + "_forceofinfection2.h5";
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
    //Uncomment to print the reproduction number
    std::cout << "The reproduction number Rc = " << sim.get_reproductionnumber_c() << "\n";

    // Uncomment to print the the values T_z1^z2
    for (int i = 0; i < (int)sim.get_T().size(); i++) {
        std::cout << "T_" << i << " = " << sim.get_T()[i] << "\n";
    }

    // Uncomment to print W_i

    std::cout << "W_i = " << sim.get_W() << "\n";

    // Uncomment to print the equilibria
    std::cout << "l_1 = " << sim.get_Equilibirum_FoI()[0] << "\n";
    std::cout << "l_2 = " << sim.get_Equilibirum_FoI()[1] << "\n";
    std::cout << "s_1 = " << sim.get_Equilibirum_compartments()[0][(int)mio::endisird::InfectionState::Susceptible]
              << "\n";
    std::cout << "s_2 = " << sim.get_Equilibirum_compartments()[1][(int)mio::endisird::InfectionState::Susceptible]
              << "\n";
    std::cout << "i_1 = " << sim.get_Equilibirum_compartments()[0][(int)mio::endisird::InfectionState::Infected]
              << "\n";
    std::cout << "i_2 = " << sim.get_Equilibirum_compartments()[1][(int)mio::endisird::InfectionState::Infected]
              << "\n";
    std::cout << "r_1 = " << sim.get_Equilibirum_compartments()[0][(int)mio::endisird::InfectionState::Recovered]
              << "\n";
    std::cout << "r_2 = " << sim.get_Equilibirum_compartments()[1][(int)mio::endisird::InfectionState::Recovered]
              << "\n";
    std::cout << "sigmaid_1 = "
              << sim.get_Equilibirum_transitions()[0][(int)mio::endisird::InfectionTransition::InfectedToDead] << "\n";
    std::cout << "sigmaid_2 = "
              << sim.get_Equilibirum_transitions()[1][(int)mio::endisird::InfectionTransition::InfectedToDead] << "\n";
    std::cout << "sigmair_1 = "
              << sim.get_Equilibirum_transitions()[0][(int)mio::endisird::InfectionTransition::InfectedToRecovered]
              << "\n";
    std::cout << "sigmair_2 = "
              << sim.get_Equilibirum_transitions()[1][(int)mio::endisird::InfectionTransition::InfectedToRecovered]
              << "\n";

    return mio::success();
}

int main()
{
    std::string result_dir = "/localdata1/trit_ha/code/memilio-1/PythonPlotsEndIDESIRD/simulation_results/";

    // Define tmax.
    ScalarType tmax = 100;

    auto result_ide = simulate_endidemodel(tmax, result_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}