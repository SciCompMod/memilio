#include "ide_endemic_secir/computed_parameters.h"
#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/infection_state.h"
#include "ide_endemic_secir/normalized_model.h"
#include "ide_endemic_secir/parameters.h"
#include "ide_endemic_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/data/analyze_result.h"
#include <Eigen/src/Core/util/Meta.h>
#include <utility>
#include <vector>

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax = 1;
    ScalarType dt   = 0.00001;

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

    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    // Uncomment for Lognorm.
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

    computed_parameters.parameters.set<mio::endisecir::NaturalBirthRate>(4e-1);
    computed_parameters.parameters.set<mio::endisecir::NaturalDeathRate>(3e-1);

    computed_parameters.set_tol_for_support_max(1e-6);

    mio::endisecir::Model model(computed_parameters);

    mio::endisecir::NormModel normmodel(computed_parameters);

    // start the simulation.
    mio::endisecir::Simulation sim(computed_parameters, model, normmodel, dt);
    sim.advance(tmax);

    //Get the compartments of model and print them.
    auto interpolated_results = mio::interpolate_simulation_result(sim.get_compartments(), dt / 2.);
    interpolated_results.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
    // Uncomment to print the compartments computed with the update scheme.
    // auto interpolated_results_update = mio::interpolate_simulation_result(sim.get_compartments_update(), dt / 2.);
    // interpolated_results_update.print_table({"US", "UE", "UC", "UI", "UH", "UU", "UR", "UD"}, 16, 8);

    //Get the commpartments of normmodel and print them.
    // auto interpolated_normresults = mio::interpolate_simulation_result(sim.get_normmodel_compartments(), dt / 2.);
    // interpolated_normresults.print_table({"s", "e", "c", "i", "h", "u", "r"}, 16, 8);

    // Uncomment to print the normalized compartments of model.
    // sim.get_normalizedcompartments().print_table({"S/N", "E/N", "C/N", "I/N", "H/N", "U/N", "R/N"}, 16, 8);

    // Uncomment to print the transitions of model.
    // sim.get_transitions().print_table({"S->E", "E->C", "C->I", "C->R", "I->H", "I->R", "H->U", "H->R", "U->D", "U->R"},
    //                                   16, 8);
    // sim.get_transitions_update().print_table(
    //     {"US->UE", "UE->UC", "UC->UI", "UC->UR", "UI->UH", "UI->UR", "UH->UU", "UH->UR", "UU->UD", "UU->UR"}, 16, 8);

    // Uncomment to print the transitions of normmodel.
    // sim.get_normmodel_transitions().print_table(
    //     {"s->e", "e->c", "c->i", "c->r", "i->h", "i->r , "h->u", "h->r", "u->d", "u->r"}, 16, 8);

    // Uncomment to print the force of infection of model.
    auto interpolated_FoI = mio::interpolate_simulation_result(sim.get_forceofinfections(), dt / 2.);
    interpolated_FoI.print_table({"FoI"}, 16, 8);
    // sim.get_forceofinfections_update().print_table({"FoIUpdate"}, 16, 8);

    // Uncomment to print the force of infection of normmodel.
    // sim.get_normmodel_forceofinfections().print_table({"norm FoI"}, 16, 8);

    //Uncomment to print the reproduction number
    std::cout << "The reproduction number Rc = " << sim.get_reproductionnumber_c() << "\n";

    // Uncomment to print the the values T_z1^z2
    //for (int i = 0; i < (int)sim.get_T().size(); i++) {
    //    std::cout << "T_" << i << " = " << sim.get_T()[i] << "\n";
    //}

    // Uncomment to print the the values V^z
    //for (int i = 0; i < (int)sim.get_V().size(); i++) {
    //    std::cout << "V_" << i << " = " << sim.get_V()[i] << "\n";
    //}

    // Uncomment to print the values W_z
    //for (int i = 0; i < (int)sim.get_W().size(); i++) {
    //    std::cout << "W_" << i << " = " << sim.get_W()[i] << "\n";
    //}

    // Uncomment to print the total population size.
    // sim.get_totalpopulations().print_table({"N"}, 16, 9);
    // sim.get_totalpopulations_update().print_table({"UN"}, 16, 9);
}