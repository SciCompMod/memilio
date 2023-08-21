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
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"
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
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(IdeOdeSecir, compareIdeOde)
{
    ScalarType t0   = 0;
    ScalarType tmax = 14.3;
    ScalarType dt   = 0.01;

    ScalarType cont_freq = 10;

    // ODE simulation

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model_ode(1);

    // Set parameters

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix
    model_ode.parameters.set<mio::osecir::Seasonality>(0.0);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    // TODO: check if this is the right way to set the contact matrix for the ODE model
    model_ode.parameters.get<mio::osecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Parameters needed to determine transition rates
    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] = 5.4;
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] = 4.0;

    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 0.4;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 0.4;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 0.4;

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
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.1;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.5;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.5;
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted=CriticalPerSevere and deathsPerSevereAdjusted = 0
    model_ode.parameters.get<mio::osecir::ICUCapacity>() = std::numeric_limits<ScalarType>::max();

    model_ode.check_constraints();

    // TODO: find out how we can change the integrator to another one from boost that doesn't use adaptive time steps
    // auto integrator = std::make_shared<mio::RKIntegratorCore>();
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // auto integrator = std::make_shared<boost::numeric::odeint::runge_kutta_cash_karp54>();
    // choose dt_min = dt_max so that we have a fixed time step and can compare to IDE
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    // set tolerance as follows so that every time step is only computed once (found out by trying)
    integrator->set_rel_tolerance(1e-1);
    integrator->set_abs_tolerance(1e-1);

    mio::TimeSeries<ScalarType> secihurd_ode = simulate(t0, tmax, dt, model_ode, integrator);

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
    ScalarType t0_ide      = 13.5;
    ScalarType Dead_before = secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 2]
                                         [(int)mio::osecir::InfectionState::Dead];

    // std::cout << "ODE numtimepoints: " << (Eigen::Index)secihurd_ode.get_num_time_points() << "\n";
    // std::cout << "Time point at t0_ide -dt: "
    //           << (Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 2 << "\n";
    // std::cout << "Dead ODE at t0_ide - dt: "
    //           << secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 2]
    //                          [(int)mio::osecir::InfectionState::Dead]
    //           << "\n";

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    mio::TimeSeries<ScalarType> init_transitions(num_transitions);
    Vec vec_init(num_transitions);
    // add dummy time point so that model initialization works
    // TODO: check if it is possible to initialize with an empty time series for init_transitions
    init_transitions.add_time_point(0, vec_init);

    // Initialize model.
    mio::isecir::Model model_ide(std::move(init_transitions), N, Dead_before);

    // Set working parameters.

    // Contact matrix; contact_matrix was already defined for ODE
    model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // To compare with the ODE model we use ExponentialDecay functions as TransitionDistributions
    // We set the funcparams so that they correspond to the above ODE model
    //TODO: We have to carefully set the default parameter or set the parameter for S to E explicitly to not have any surprises
    // when calculating the calc_time when computing compartments
    // maybe we should check this in check constraints or set this support_max explicitly to 0 so that the user
    // doesn't have to think of this
    // Why do we not need the distribution from S to E? Compartment S and the flow from S to E is calculated in another way and not via
    // the compute_compartments/compute_flows Funktion.
    mio::ExponentialDecay expdecay(10.0);
    mio::StateAgeFunctionWrapper delaydistribution(expdecay);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    // TransitionDistribution from S to E is never used, use the same one as in TransitionDistribution from E to C
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_parameter(
        (2 * model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] -
         model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]));
    // see definition of rate_E in model.h of ODE; set funcparam to 1/rate_E
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_parameter(
        (2 * model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] -
         model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]));
    std::cout << "rate_E: "
              << 1 / (2 * model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] -
                      model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0])
              << "\n";
    // see definition of rate_E in model.h of ODE; set funcparam to 1/rate_INS
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms].set_parameter(
        2 * (model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] -
             model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]));
    std::cout << "rate_INS: "
              << 1 / (2 * (model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] -
                           model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]))
              << "\n";
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_parameter(
        2 * (model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] -
             model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]));
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_parameter(
        1 / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_parameter(
        1 / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical].set_parameter(
        1 / model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_parameter(
        1 / model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_parameter(
        1 / model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_parameter(
        1 / model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);

    model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set probabilities that determine proportion between compartments
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    // // Set contact matrix
    // mio::ContactMatrixGroup contact_matrix_ide = mio::ContactMatrixGroup(1, 1);
    // contact_matrix_ide[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    // model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix_ide);

    // Set further parameters
    mio::ConstantFunction constfunc_proboncontact(
        model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
    mio::StateAgeFunctionWrapper proboncontact(constfunc_proboncontact);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(proboncontact);

    mio::ConstantFunction constfunc_reltransnosy(
        model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
    mio::StateAgeFunctionWrapper reltransnosy(constfunc_reltransnosy);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(reltransnosy);

    mio::ConstantFunction constfunc_riskofinf(
        model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
    mio::StateAgeFunctionWrapper riskofinf(constfunc_riskofinf);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(riskofinf);

    // Set t0 for IDE simulation
    // ScalarType t0_ide = 6 * model_ide.get_global_support_max(dt);
    std::cout << "t0_ide: " << t0_ide << "\n";

    // Compute initial flows from results of ODE simulation
    model_ide.compute_initial_flows_from_compartments2(secihurd_ode, t0_ide, dt);

    model_ide.check_constraints(dt);

    std::cout << model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::Susceptible)]
              << "\n";

    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::Susceptible)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::Susceptible)];
    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::Exposed)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::Exposed)];
    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::InfectedNoSymptoms)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::InfectedNoSymptoms)];
    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::InfectedSymptoms)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::InfectedSymptoms)];
    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::InfectedSevere)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::InfectedSevere)];
    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::InfectedCritical)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::InfectedCritical)];
    model_ide.m_populations.get_last_value()[Eigen::Index(mio::isecir::InfectionState::Recovered)] =
        secihurd_ode[(t0_ide - t0) / dt][Eigen::Index(mio::isecir::InfectionState::Recovered)];

    // Carry out simulation
    std::cout << "Simulating now \n";
    mio::isecir::Simulation sim(model_ide, t0_ide, dt);
    sim.advance(tmax);
    // sim.print_transitions();
    sim.print_compartments();
    // sim.print_transitions();
    // mio::TimeSeries<ScalarType> secihurd_ide = sim.get_result();
}