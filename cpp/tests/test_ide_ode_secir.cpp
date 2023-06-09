/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
// #include "Eigen/src/Core/util/Meta.h"
#include "matchers.h"
// #include "load_test_data.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"
#include "memilio/math/adapt_rk.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameters.h"
#include <distributions_helpers.h>
#include <gtest/gtest.h>

#include "boost/fusion/functional/invocation/invoke.hpp"
#include "load_test_data.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <gtest/gtest.h>

// define function that computes flows needed for initalization of IDE for a given result/compartments of the ODE model
// we assume that the ODE simulation starts at t0=0
void compute_initial_flows_from_compartments(mio::TimeSeries<ScalarType> secihurd_ode, mio::isecir::Model model,
                                             ScalarType t0_ide, ScalarType dt)
{
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // get (global) max_support to determine how many flows in the past we have to compute
    ScalarType global_max_support = model.get_global_max_support(dt);
    std::cout << "Global max_support: " << global_max_support << "\n";
    Eigen::Index global_max_support_index = std::ceil(global_max_support / dt);

    // remove time point
    model.m_transitions.remove_last_time_point();

    ScalarType t0_ide_index = std::ceil(t0_ide / dt);
    unused(secihurd_ode);
    // flow from S to E for -6*global_max_support, ..., 0 (directly from compartments)
    // add time points to init_transitions here
    for (int i = t0_ide_index - 6 * global_max_support_index + 1; i <= t0_ide_index; i++) {
        model.m_transitions.add_time_point(i, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        std::cout << "i: " << i << "\n";
        model.m_transitions.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] =
            secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::Susceptible)] -
            secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::Susceptible)];
    }

    // then use compute_flow function to compute following flows

    Eigen::Index start_shift = t0_ide_index - 6 * global_max_support_index;
    std::cout << "first timepoint: " << start_shift << "\n";

    // flow from E to C for -5*global_max_support, ..., 0
    for (int i = t0_ide_index - 5 * global_max_support_index + 1; i <= t0_ide_index; i++) {
        model.compute_flow(1, 0, dt, true, i - start_shift);
    }

    // flow from C to I and C to R for -4*global_max_support, ..., 0
    for (int i = t0_ide_index - 4 * global_max_support_index + 1; i <= t0_ide_index; i++) {
        // C to I
        model.compute_flow(2, 1, dt, true, i - start_shift);
        // C to R
        model.compute_flow(3, 1, dt, true, i - start_shift);
    }

    // flow from I to H and I to R for -3*global_max_support, ..., 0
    for (int i = t0_ide_index - 3 * global_max_support_index + 1; i <= t0_ide_index; i++) {
        // I to H
        model.compute_flow(4, 2, dt, true, i - start_shift);
        // I to R
        model.compute_flow(5, 2, dt, true, i - start_shift);
    }

    // flow from H to U and H to R for -2*global_max_support, ..., 0
    for (int i = t0_ide_index - 2 * global_max_support_index + 1; i <= t0_ide_index; i++) {
        // H to U
        model.compute_flow(6, 4, dt, true, i - start_shift);
        // H to U
        model.compute_flow(7, 4, dt, true, i - start_shift);
    }

    // flow from U to D and U to R for -1*global_max_support, ..., 0
    for (int i = t0_ide_index - 1 * global_max_support_index + 1; i <= t0_ide_index; i++) {
        // U to D
        model.compute_flow(8, 6, dt, true, i - start_shift);
        // U to R
        model.compute_flow(9, 6, dt, true, i - start_shift);
    }

    std::cout << "numtimepoints: " << model.m_transitions.get_num_time_points() << "\n";

    mio::isecir::Simulation sim(model, t0_ide, dt);
    sim.print_transitions();

    // return TimeSeries from -1*global_max_support, ..., 0 containing necessary flows for IDE simulation
}

TEST(IdeOdeSecir, compareIdeOde)
{
    ScalarType t0   = 0;
    ScalarType tmax = 25;
    ScalarType dt   = 1;

    ScalarType cont_freq = 10;

    // ODE simulation

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model_ode(1);

    //Set parameters

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix
    model_ode.parameters.set<mio::osecir::Seasonality>(0.0);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    // Parameters needed to determine transition rates
    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] = 5.2;
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] = 4.0;

    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 2.0;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 2.0;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 2.0;

    // Set initial values for compartments
    model_ode.populations.set_total(nb_total_t0);
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    nb_total_t0);

    // Set probabilities that determine proportion between compartments
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.5;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 0.5;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 0.5;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.5;

    // Further model parameters
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.5;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.5;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.5;
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<double>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted=CriticalPerSevere and deathsPerSevereAdjusted = 0
    model_ode.parameters.get<mio::osecir::ICUCapacity>() = std::numeric_limits<double>::max();

    model_ode.check_constraints();

    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    // choose dt_min = dt_max so that we have a fixed time step and can compare to IDE
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    // set tolerance as follows so that every time step is only computed once (found out by trying)
    integrator->set_rel_tolerance(1e-2);
    integrator->set_abs_tolerance(1e-2);

    mio::TimeSeries<double> secihurd_ode = simulate(t0, tmax, dt, model_ode, integrator);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        char vars[] = {'S', 'E', 'C', 'I', 'H', 'U', 'R', 'D'};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::osecir::InfectionState::Count; k++) {
            printf(" %c", vars[k]);
        }
        auto num_points = static_cast<size_t>(secihurd_ode.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", secihurd_ode.get_time(i));
            Eigen::VectorXd res_j = secihurd_ode.get_value(i);
            for (size_t j = 0; j < (size_t)mio::osecir::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }
        std::cout << "\n";
        Eigen::VectorXd res_j = secihurd_ode.get_last_value();
        printf("number total: %f",
               res_j[0] + res_j[1] + res_j[2] + res_j[3] + res_j[4] + res_j[5] + res_j[6] + res_j[7]);
    }

    std::cout << "\n";

    // IDE simulation

    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType N           = nb_total_t0;
    ScalarType Dead_before = 10;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    mio::TimeSeries<ScalarType> init_transitions(num_transitions);
    Vec vec_init(num_transitions);
    init_transitions.add_time_point(0, vec_init);

    // Initialize model.
    mio::isecir::Model model_ide(std::move(init_transitions), N, Dead_before);

    ScalarType t0_ide = 6 * model_ide.get_global_max_support(dt) + 3;
    std::cout << "t0_ide: " << t0_ide << "\n";
    compute_initial_flows_from_compartments(secihurd_ode, model_ide, t0_ide, dt);

    // Set working parameters.
    // To compare woth the ODE model we use ExponentialDecay functions as TransitionDistributions
    // We set the funcparams so that they correspond to the above ODE model
    mio::isecir::ExponentialDecay expdecay(4);
    mio::isecir::StateAgeFunctionWrapper delaydistribution;
    delaydistribution.set_state_age_function(expdecay);
    std::vector<mio::isecir::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set probabilities that determine proportion between compartments
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix_ide = mio::ContactMatrixGroup(1, 1);
    contact_matrix_ide[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix_ide);

    mio::isecir::StateAgeFunctionWrapper prob;
    mio::isecir::ConstantFunction constfunc(0.5);
    prob.set_state_age_function(constfunc);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    model_ide.check_constraints(dt);

    // Carry out simulation.
    // mio::isecir::Simulation sim(model_ide, t0_ide - 1, dt);
    // sim.advance(tmax);
    // mio::TimeSeries<ScalarType> secihurd_ide = sim.get_result();
}