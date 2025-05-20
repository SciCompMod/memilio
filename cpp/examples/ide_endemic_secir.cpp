#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/infection_state.h"
#include "ide_endemic_secir/parameters.h"
#include "ide_endemic_secir/simulation.h"
#include "ide_endemic_secir/infection_state.h"
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

    ScalarType tmax = 10;
    ScalarType dt   = 0.1;

    int num_states      = static_cast<int>(mio::endisecir::InfectionState::Count);
    int num_transitions = static_cast<int>(mio::endisecir::InfectionTransition::Count);

    // Create TimeSeries with num_states elements where states needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_states);

    Vec vec_init(num_states);

    vec_init[static_cast<int>(mio::endisecir::InfectionState::Susceptible)]        = 9950.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Exposed)]            = 25.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedNoSymptoms)] = 15.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedSymptoms)]   = 8.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedSevere)]     = 1.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::InfectedCritical)]   = 1.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Recovered)]          = 5.;
    vec_init[static_cast<int>(mio::endisecir::InfectionState::Dead)]               = 0.;

    init.add_time_point(0, vec_init);

    mio::endisecir::Model model(std::move(init));

    //Set working parameters

    // mio::ExponentialSurvivalFunction exp(2.0);
    // mio::StateAgeFunctionWrapper delaydistribution(exp);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    mio::SmootherCosine smoothcos(4.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistribution(num_transitions, delaydistribution);

    model.parameters.get<mio::endisecir::TransitionDistributions>() = vec_delaydistribution;

    std::vector<ScalarType> vec_prob(num_transitions, 0.5);
    // The following probabilities must be 1, as there is no other way to go.
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::endisecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.get<mio::endisecir::TransitionProbabilities>()                          = vec_prob;

    mio::ContactMatrixGroup contact_matrix                  = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
    model.parameters.get<mio::endisecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction constant(0.1);
    mio::StateAgeFunctionWrapper constant_prob(constant);

    model.parameters.get<mio::endisecir::TransmissionProbabilityOnContact>() = constant_prob;

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper exponential_prob(exponential);

    model.parameters.get<mio::endisecir::RelativeTransmissionNoSymptoms>() = exponential_prob;
    model.parameters.get<mio::endisecir::RiskOfInfectionFromSymptomatic>() = exponential_prob;

    model.parameters.set<mio::endisecir::NaturalBirthRate>(5e-3);
    model.parameters.set<mio::endisecir::NaturalDeathRate>(3e-3);

    // start the simulation.
    mio::endisecir::Simulation sim(model, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_compartments(), dt / 2.);

    interpolated_results.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);

    // Uncomment to print the reproduction number
    std::cout << "The reproduction number Rc = " << sim.get_reproductionnumber_c() << "\n";

    // Uncomment to print the transitions.
    // sim.get_transitions().print_table(
    //     {"S->E 1", "E->C 1", "C->I 1", "C->R 1", "I->H 1", "I->R 1", "H->U 1", "H->R 1", "U->D 1", "U->R 1"}, 16, 8);

    // Uncomment to print the normalized compartments.
    // sim.get_normalizedcompartments().print_table({"s", "e", "c", "i", "h", "u", "r", "d "}, 16, 8);

    // Uncomment to print the total population size.
    sim.get_totalpopulations().print_table({"N"}, 16, 9);

    // Uncomment to print the force of infection.
    sim.get_forceofinfections().print_table({"FoI"}, 16, 8);

    // std::vector<ScalarType> equi = sim.get_equilibriumcompartments();
    // std::cout << "Equilibrium normalized compartments: \n";
    // std::cout << "foi* " << sim.get_equilibrium_forceofinfection() << "\n";
    // std::cout << "s* " << equi[(int)mio::endisecir::InfectionState::Susceptible] << "\n";
    // std::cout << "e* " << equi[(int)mio::endisecir::InfectionState::Exposed] << "\n";
    // std::cout << "c* " << equi[(int)mio::endisecir::InfectionState::InfectedNoSymptoms] << "\n";
    // std::cout << "i* " << equi[(int)mio::endisecir::InfectionState::InfectedSymptoms] << "\n";
    // std::cout << "h* " << equi[(int)mio::endisecir::InfectionState::InfectedSevere] << "\n";
    // std::cout << "u* " << equi[(int)mio::endisecir::InfectionState::InfectedCritical] << "\n";
    // std::cout << "r* " << equi[(int)mio::endisecir::InfectionState::Recovered] << "\n";
}