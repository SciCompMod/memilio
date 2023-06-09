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

    // // Set working parameters.
    // // To compare woth the ODE model we use ExponentialDecay functions as TransitionDistributions
    // // We set the funcparams so that they correspond to the above ODE model
    // mio::isecir::ExponentialDecay expdecay(4);
    // mio::isecir::StateAgeFunctionWrapper delaydistribution;
    // delaydistribution.set_state_age_function(expdecay);
    // std::vector<mio::isecir::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    // model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

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

    ScalarType t0_ide = 6 * model_ide.get_global_max_support(dt) + 3;
    std::cout << "t0_ide: " << t0_ide << "\n";

    model_ide.compute_initial_flows_from_compartments(secihurd_ode, t0_ide, dt);

    std::cout << "m_transitions num_time_points: " << model_ide.m_transitions.get_num_time_points() << "\n";
    std::cout << "global_max_supp_timesteps: " << (Eigen::Index)std::ceil(model_ide.get_global_max_support(dt) / dt)
              << "\n";

    model_ide.check_constraints(dt);

    // Carry out simulation.
    std::cout << "Simulating now \n";
    mio::isecir::Simulation sim(model_ide, t0_ide, dt);

    sim.advance(tmax);
    sim.print_compartments();
    mio::TimeSeries<ScalarType> secihurd_ide = sim.get_result();
}