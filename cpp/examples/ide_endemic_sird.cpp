#include "ide_endemic_sird/computed_parameters.h"
#include "ide_endemic_sird/model.h"
#include "ide_endemic_sird/infection_state.h"
#include "ide_endemic_sird/normalized_model.h"
#include "ide_endemic_sird/parameters.h"
#include "ide_endemic_sird/simulation.h"
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

    ScalarType tmax = 50;
    ScalarType dt   = 1.0;

    int num_states      = static_cast<int>(mio::endisird::InfectionState::Count);
    int num_transitions = static_cast<int>(mio::endisird::InfectionTransition::Count);

    // Create TimeSeries with num_states elements where states needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_states);

    Vec vec_init(num_states);

    vec_init[static_cast<int>(mio::endisird::InfectionState::Susceptible)] = 100000.;
    vec_init[static_cast<int>(mio::endisird::InfectionState::Infected)]    = 40.;
    vec_init[static_cast<int>(mio::endisird::InfectionState::Recovered)]   = 0.;
    vec_init[static_cast<int>(mio::endisird::InfectionState::Dead)]        = 0.;

    init.add_time_point(0, vec_init);

    mio::endisird::CompParameters computed_parameters(std::move(init));

    //Set working parameters

    // mio::ExponentialSurvivalFunction exp(3.0);
    // mio::StateAgeFunctionWrapper delaydistribution(exp);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    // mio::SmootherCosine smoothcos(2.0);
    // mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    mio::ConstantFunction initialfunc(0);
    mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistributioninit);
    // InfectedToDead
    mio::SmootherCosine survivalExposedToInfectedNoSymptoms(6.0);
    vec_delaydistribution[(int)mio::endisird::InfectionTransition::InfectedToDead].set_state_age_function(
        survivalExposedToInfectedNoSymptoms);
    // InfectedToRecovered
    mio::SmootherCosine survivalInfectedNoSymptomsToInfectedSymptoms(8.0);
    vec_delaydistribution[(int)mio::endisird::InfectionTransition::InfectedToRecovered].set_state_age_function(
        survivalInfectedNoSymptomsToInfectedSymptoms);

    // mio::ConstantFunction initialfunc(0);
    // mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistributioninit);
    // // InfectedToDead
    // mio::LognormSurvivalFunction survivalExposedToInfectedNoSymptoms(0.8, 0, 4.2);
    // vec_delaydistribution[(int)mio::endisird::InfectionTransition::InfectedToDead].set_state_age_function(
    //     survivalExposedToInfectedNoSymptoms);
    // // InfectedToRecovered
    // mio::LognormSurvivalFunction survivalInfectedNoSymptomsToInfectedSymptoms(0.2, 0, 0.8);
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

    computed_parameters.parameters.set<mio::endisird::NaturalBirthRate>(4e-3);
    computed_parameters.parameters.set<mio::endisird::NaturalDeathRate>(3e-3);

    computed_parameters.set_tol_for_support_max(1e-6);

    mio::endisird::Model model(computed_parameters);

    mio::endisird::NormModel normmodel(computed_parameters);

    // start the simulation.
    mio::endisird::Simulation sim(computed_parameters, model, normmodel, dt);
    sim.advance(tmax);

    //Get the compartments of model and print them.
    auto interpolated_results = mio::interpolate_simulation_result(sim.get_compartments(), dt / 2.);
    interpolated_results.print_table({"S", "I", "R", "D "}, 16, 8);

    //Get the commpartments of normmodel and print them.
    // auto interpolated_normresults = mio::interpolate_simulation_result(sim.get_normmodel_compartments(), dt / 2.);
    // interpolated_normresults.print_table({"s", "i",  "r"}, 16, 8);

    // Uncomment to print the normalized compartments of model.
    // sim.get_normalizedcompartments().print_table({"S/N", "I/N",  "R/N"}, 16, 8);

    // Uncomment to print the transitions of model.
    // sim.get_transitions().print_table({"S->I", "I->D", "I->R"},
    //                                   16, 8);

    // Uncomment to print the transitions of normmodel.
    // sim.get_normmodel_transitions().print_table(
    //     {"s->i", "i->d", "i->r }, 16, 8);

    // Uncomment to print the force of infection of model.
    auto interpolated_FoI = mio::interpolate_simulation_result(sim.get_forceofinfections(), dt / 2.);
    interpolated_FoI.print_table({"FoI"}, 16, 8);

    // Uncomment to print the force of infection of normmodel.
    // sim.get_normmodel_forceofinfections().print_table({"norm FoI"}, 16, 8);

    //Uncomment to print the reproduction number
    std::cout << "The reproduction number Rc = " << sim.get_reproductionnumber_c() << "\n";

    // Uncomment to print the the values T_z1^z2
    //for (int i = 0; i < (int)sim.get_T().size(); i++) {
    //    std::cout << "T_" << i << " = " << sim.get_T()[i] << "\n";
    //}

    // Uncomment to print the values W_z
    //for (int i = 0; i < (int)sim.get_W().size(); i++) {
    //    std::cout << "W_" << i << " = " << sim.get_W()[i] << "\n";
    //}

    // Uncomment to print the total population size.
    // sim.get_totalpopulations().print_table({"N"}, 16, 9);
}