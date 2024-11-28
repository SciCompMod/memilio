#include "load_test_data.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <gtest/gtest.h>
#include <vector>

TEST(TestIdeAgeres, compareWithPreviousRun)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    size_t num_agegroups = 3;
    ScalarType tmax      = 5.0;
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 5000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 6.);
    ScalarType dt = 1.;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions * num_agegroups elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions * num_agegroups);

    // Define transitions that will be used for initialization.
    Vec vec_init = Vec::Constant(num_transitions * num_agegroups, 1.);
    // First AgeGroup.
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 20.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 3.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 10.0;
    // Second AgeGroup.
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 20.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 10.0;
    // Third AgeGroup.
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
    // Add initial time point to TimeSeries.
    init.add_time_point(-10., vec_init);
    // Add further time points until t0.
    while (init.get_last_time() < -dt / 2.) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths, num_agegroups);

    // Set working parameters.

    // First AgeGroup for TransitionDistributions.
    mio::SmootherCosine smoothcos1(2.0);
    mio::StateAgeFunctionWrapper delaydistribution1(smoothcos1);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib1(num_transitions, delaydistribution1);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(3.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(0)] = vec_delaydistrib1;

    // Second AgeGroup for TransitionDistributions.
    mio::SmootherCosine smoothcos2(3.0);
    mio::StateAgeFunctionWrapper delaydistribution2(smoothcos2);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib2(num_transitions, delaydistribution2);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(2.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(1)] = vec_delaydistrib2;

    // Third AgeGroup for TransitionDistributions.
    mio::SmootherCosine smoothcos3(2.5);
    mio::StateAgeFunctionWrapper delaydistribution3(smoothcos3);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib3(num_transitions, delaydistribution3);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(4.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(2)] = vec_delaydistrib3;

    //TransitionProbabilities.
    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        std::vector<ScalarType> vec_prob(num_transitions, 0.5);
        // The following probabilities must be 1, as there is no other way to go.
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        model.parameters.get<mio::isecir::TransitionProbabilities>()[group]                   = vec_prob;
    }

    // Contact matrix.
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, num_agegroups);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper prob(exponential);
    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = prob;
        model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = prob;
        model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = prob;
    }
    model.parameters.set<mio::isecir::Seasonality>(0.1);
    // Start the simulation on the 40th day of a year (i.e. in February).
    model.parameters.set<mio::isecir::StartDay>(40);

    model.check_constraints(dt);

    //Compare compartments with previous run.
    mio::TimeSeries<ScalarType> compartments = simulate(tmax, dt, model);
    auto compare_compartments                = load_test_data_csv<ScalarType>("ide-secir-ageres-compare.csv");

    ASSERT_EQ(compare_compartments.size(), static_cast<size_t>(compartments.get_num_time_points()));
    for (size_t i = 0; i < compare_compartments.size(); i++) {
        ASSERT_EQ(compare_compartments[i].size(), static_cast<size_t>(compartments.get_num_elements()) + 1)
            << "at row " << i;
        ASSERT_NEAR(compartments.get_time(i), compare_compartments[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare_compartments[i].size(); j++) {
            ASSERT_NEAR(compartments.get_value(i)[j - 1], compare_compartments[i][j], 1e-7) << " at row " << i;
        }
    }
    //Compare transitions with previous run.

    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);

    auto transitions         = sim.get_transitions();
    auto compare_transitions = load_test_data_csv<ScalarType>("ide-secir-ageres-transitions-compare.csv");

    size_t iter_0 = 0;
    while (transitions.get_time(iter_0) < compare_transitions[0][0]) {
        iter_0++;
    }

    for (size_t i = 0; i < compare_transitions.size(); i++) {
        ASSERT_EQ(compare_transitions[i].size(), static_cast<size_t>(transitions.get_num_elements()) + 1)
            << "at row " << i;
        ASSERT_NEAR(transitions.get_time(i + iter_0), compare_transitions[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare_transitions[i].size(); j++) {
            ASSERT_NEAR(transitions.get_value(i + iter_0)[j - 1], compare_transitions[i][j], 1e-7) << " at row " << i;
        }
    }
}