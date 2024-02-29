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
// #include "load_test_data.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"
#include "memilio/math/adapt_rk.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameters.h"

#include "boost/fusion/functional/invocation/invoke.hpp"
#include "ide_secir/infection_state.h"
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "ide_secir/initialize_from_ode.h"
#include "memilio/io/result_io.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <string>

int main()
{

    // Here we decide what exactly we want to in the example below
    bool print_to_terminal = false;
    bool save_result       = true;
    bool ide_simulation    = true;
    int dt_exponent        = 3;
    // We use setting 2 as baseline, changes for other settings are in respective if statements
    int setting = 16;

    // General set up.
    ScalarType t0   = 0;
    ScalarType tmax = 70.00;
    ScalarType dt   = pow(10, -dt_exponent);

    /**********************************
    *         ODE simulation          *
    **********************************/

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 20, nb_car_t0 = 20, nb_inf_t0 = 3, nb_hosp_t0 = 1, nb_icu_t0 = 1,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    if (setting == 10 || setting == 12 || setting == 13 || setting == 14) {
        nb_rec_t0 = 0.;
    }

    mio::osecir::Model model_ode(1);

    // Set parameters
    ScalarType cont_freq = 1.0;

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix
    model_ode.parameters.set<mio::osecir::Seasonality>(0.0);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_ode.parameters.get<mio::osecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Parameters needed to determine transition rates
    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] = 2.6;
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] = 2.0;

    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 0.3;

    // Set initial values for compartments
    model_ode.populations.set_total(nb_total_t0);
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = nb_exp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = nb_car_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = nb_inf_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = nb_hosp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = nb_icu_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = nb_rec_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = nb_dead_t0;
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    nb_total_t0);

    // Set probabilities that determine proportion between compartments
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.5;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 0.5;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 0.5;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.5;

    if (setting == 4) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.5;
    }

    if (setting == 5 || setting == 12 || setting == 13 || setting == 15) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 1.;
    }

    if (setting == 6 || setting == 11 || setting == 14) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.;
    }

    if (setting == 7) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 1.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 0.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 0.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.;
    }

    if (setting == 8) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0] = 0.;
    }

    // Further model parameters
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 1.0;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 1.0;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 1.0;
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0
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
               res_j[0] + res_j[1] + res_j[2] + res_j[4] + res_j[6] + res_j[7] + res_j[8] + res_j[9]);
    }

    std::cout << "\n";

    if (save_result) {
        auto save_result_status_ode = mio::save_result({secihurd_ode}, {0}, 1,
                                                       "../../results/result_ode_dt=1e-" + std::to_string(dt_exponent) +
                                                           "_setting" + std::to_string(setting) + ".h5");
    }

    /**********************************
    *         IDE simulation          *
    **********************************/

    if (ide_simulation) {
        using Vec = mio::TimeSeries<ScalarType>::Vector;

        ScalarType N = nb_total_t0;
        // TODO: Set this automatically or check if this is possible (wrt to global_max_support) with the given ODE simulation
        ScalarType t0_ide = 35.0;
        // Get number of dead individuals at time -dt from ODE model
        ScalarType deaths = secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                                        [(int)mio::osecir::InfectionState::Dead];

        int num_transitions = (int)mio::isecir::InfectionTransition::Count;

        mio::TimeSeries<ScalarType> init_transitions(num_transitions);
        Vec vec_init(num_transitions);
        // Add dummy time point so that model initialization works
        // Attention: here we need to initilaize with a time series that last time point is t0_ide so that m_populations is set correctly
        // TODO: check if it is possible to initialize with an empty time series for init_transitions
        // TODO: check if there is an easier/cleaner way to initialize
        init_transitions.add_time_point(t0_ide, vec_init);

        ScalarType total_infections = 0.;

        if (setting == 11 || setting == 13 || setting == 14) {
            total_infections = secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                                           [(int)mio::osecir::InfectionState::InfectedSymptoms] +
                               secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                                           [(int)mio::osecir::InfectionState::InfectedSevere] +
                               secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                                           [(int)mio::osecir::InfectionState::InfectedCritical] +
                               secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                                           [(int)mio::osecir::InfectionState::Recovered] +
                               secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                                           [(int)mio::osecir::InfectionState::Dead];
        }

        bool need_flow_init = true;

        // Initialize model.
        mio::isecir::Parameters parameters;
        mio::isecir::Model model_ide(std::move(init_transitions), N, deaths, total_infections, parameters,
                                     need_flow_init);

        // Set working parameters.

        // Contact matrix; contact_matrix was already defined for ODE
        model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        // To compare with the ODE model we use ExponentialDecay functions as TransitionDistributions
        // We set the parameters so that they correspond to the above ODE model
        mio::ExponentialDecay expdecay(10.0);
        mio::StateAgeFunctionWrapper delaydistribution(expdecay);
        std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
        // ExposedToInfectedNoSymptoms
        // see definition of rate_E in model.h of ODE; set parameter to rate_E
        ScalarType rate_E = 1 / (2 * model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] -
                                 model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]);
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_parameter(rate_E);
        // InfectedNoSymptomsToInfectedSymptoms
        // see definition of rate_INS in model.h of ODE; set parameter to rate_INS
        ScalarType rate_INS = 1 / (2 * (model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] -
                                        model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]));
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms].set_parameter(
            rate_INS);
        // InfectedNoSymptomsToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_parameter(rate_INS);
        // InfectedSymptomsToInfectedSevere
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
        // InfectedSymptomsToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
        // InfectedSevereToInfectedCritical
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical].set_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
        // InfectedSevereToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
        // InfectedCriticalToDead]
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);
        // InfectedCriticalToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);

        model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

        // Set probabilities that determine proportion between compartments
        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0];
        model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

        // Set further parameters
        mio::ConstantFunction constfunc_proboncontact(
            model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper proboncontact(constfunc_proboncontact);
        model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(proboncontact);

        mio::ConstantFunction constfunc_reltransnosympt(
            model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper reltransnosympt(constfunc_reltransnosympt);
        model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(reltransnosympt);

        mio::ConstantFunction constfunc_riskofinf(
            model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper riskofinf(constfunc_riskofinf);
        model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(riskofinf);

        // Compute initial flows from results of ODE simulation
        compute_initial_flows_from_ode_compartments(model_ode, model_ide, secihurd_ode, t0_ide, dt);

        model_ide.check_constraints(dt);

        model_ide.set_populations_before_simulation();

        if (setting == 15 || setting == 16) {
            model_ide.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] =
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt - 1]
                            [(int)mio::osecir::InfectionState::Recovered];
        }

        // Carry out simulation
        std::cout << "Simulating now \n";
        mio::isecir::Simulation sim(model_ide, t0_ide, dt);
        sim.advance(tmax);

        mio::TimeSeries<ScalarType> secihurd_ide       = sim.get_result();
        mio::TimeSeries<ScalarType> secihurd_ide_flows = sim.get_transitions();

        if (print_to_terminal) {
            secihurd_ide.print_table();
        }

        std::cout << "Initialization method: " << sim.get_model().get_initialization_method() << "\n";

        std::cout << "Compartments at last time step of ODE:\n";
        std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
        for (Eigen::Index j = 0; j < secihurd_ode.get_num_elements(); ++j) {
            std::cout << "  |  " << std::fixed << std::setprecision(8) << secihurd_ode.get_last_value()[j];
        }
        std::cout << "\n" << std::endl;

        std::cout << "Compartments at last time step of IDE:\n";
        std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
        for (Eigen::Index j = 0; j < secihurd_ide.get_num_elements(); ++j) {
            std::cout << "  |  " << std::fixed << std::setprecision(8) << secihurd_ide.get_last_value()[j];
        }
        std::cout << "\n" << std::endl;

        if (save_result) {

            auto save_result_status_ide =
                mio::save_result({secihurd_ide}, {0}, 1,
                                 "../../results/result_ide_dt=1e-" + std::to_string(dt_exponent) + "_setting" +
                                     std::to_string(setting) + ".h5");
            // auto save_result_status_ide_flows =
            //     mio::save_result({secihurd_ide_flows}, {0}, 1,
            //                      "../../results/result_ide_flows_dt=1e-" + std::to_string(dt_exponent) + "_setting" +
            //                          std::to_string(setting) + ".h5");
        }
    }
}
